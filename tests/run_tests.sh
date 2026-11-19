#!/usr/bin/env bash
# bandLD tests: exact reference from the generator's dosage numerators (needs python3, numpy, scipy).
# usage: tests/run_tests.sh BIN        (WORKDIR=... to choose the scratch directory)
set -uo pipefail
BIN=$(readlink -f "$1"); HERE=$(cd "$(dirname "$0")" && pwd)
W=${WORKDIR:-/tmp/bandLD-tests}; mkdir -p "$W"; cd "$W"
GEN="python3 $HERE/make_bgen.py"; CHK="python3 $HERE/check.py"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "PASS  $*"; }
bad() { fail=$((fail+1)); echo "FAIL  $*"; }
gen() { local n=$1 N=$2 M=$3 s=$4; shift 4; [ -f $n.bgen ] || env "$@" $GEN $W/$n $N $M $s >/dev/null || echo "GENERATOR FAILED $n"; }

# data: positions are spread over 50*M bp, so a 4,000 bp window keeps roughly a quarter of the pairs
gen imp 1001 300 21 PROFILE=imputed FLIP=0.4
gen dense 600 300 22 PROFILE=dense
gen hard 500 257 23 PROFILE=hardcall
gen tiny 65 33 24 PROFILE=dense

for d in imp dense hard tiny; do
  for win in 4000 100000000; do
    for th in 1 4; do
      "$BIN" build --exact --bgen $d.bgen --out $d.w$win.t$th --window-bp $win --threads $th >$d.build.log 2>&1 || { bad "$d build window=$win threads=$th: $(tail -2 $d.build.log)"; continue; }
      "$BIN" query --band $d.w$win.t$th.bld --out $d.w$win.t$th.ld --min-bp 0 --max-bp 2000000000 --precision 17 --threads $th >$d.query.log 2>&1 || { bad "$d query: $(tail -2 $d.query.log)"; continue; }
      if out=$($CHK $W/$d $d.w$win.t$th.ld $win 0); then ok "$d exact window=$win threads=$th: ${out#*: }"; else echo "$out"; bad "$d exact window=$win threads=$th"; fi
    done
    cmp -s $d.w$win.t1.bld $d.w$win.t4.bld && ok "$d window=$win: 1 and 4 threads give the same file" || bad "$d window=$win: thread counts differ"
    for k in avx2 scalar; do
      BANDLD_KERNEL=$k "$BIN" build --exact --bgen $d.bgen --out $d.k --window-bp $win >/dev/null 2>&1
      cmp -s $d.w$win.t1.bld $d.k.bld && ok "$d window=$win kernel=$k: same file" || bad "$d window=$win kernel=$k differs"
    done
  done
done

# rounded mode
"$BIN" build --bgen imp.bgen --out imp.r15 --window-bp 4000 >/dev/null 2>&1          # default: 16-bit precision
"$BIN" query --band imp.r15.bld --out imp.r15.ld --min-bp 0 --max-bp 2000000000 --precision 17 >/dev/null 2>&1
if out=$($CHK $W/imp imp.r15.ld 4000 15); then ok "default precision (16-bit): ${out#*: }"; else echo "$out"; bad "rounded mode"; fi
[ $(stat -c %s imp.r15.bld) -lt $(stat -c %s imp.w4000.t1.bld) ] && ok "rounded file is smaller ($(stat -c %s imp.r15.bld) < $(stat -c %s imp.w4000.t1.bld) bytes)" || bad "rounded file not smaller"

# region query = the same block of the full matrix; gzip output
P1=$(sed -n 60p imp.z | cut -d' ' -f3); P2=$(sed -n 200p imp.z | cut -d' ' -f3)
"$BIN" query --band imp.w4000.t1.bld --out imp.reg.ld.gz --chr 01 --min-bp $P1 --max-bp $P2 --precision 17 >/dev/null 2>&1
if out=$($CHK $W/imp imp.reg.ld.gz 4000 0); then ok "region query, gzip: ${out#*: }"; else echo "$out"; bad "region query"; fi
EXP=$(awk -v a=$P1 -v b=$P2 '!/^#/ && $2 >= a && $2 <= b' imp.w4000.t1.snps.tsv | wc -l)
[ "$(grep -vc '^#' imp.reg.ld.gz.snps.tsv)" = "$EXP" ] && ok "region has the expected $EXP variants" || bad "region size: $(grep -vc '^#' imp.reg.ld.gz.snps.tsv) instead of $EXP"

# SNP list, plain IDs in a shuffled order
python3 - "$W" <<'E'
import random, sys
w = sys.argv[1]; z = [l.split() for l in open(f"{w}/imp.z").read().strip().split("\n")[1:]]
stored = {l.split("\t")[2] for l in open(f"{w}/imp.w4000.t1.snps.tsv") if not l.startswith("#")}
z = [r for r in z if r[0] in stored]          # variants the build skipped (monomorphic) cannot be queried
random.seed(5); pick = random.sample(z, 40)
open(f"{w}/list.ids", "w").write("".join(r[0] + "\n" for r in pick))
# table with alleles, every second SNP listed with swapped alleles
t = ["rsid chromosome position allele1 allele2"]
for k, r in enumerate(pick): t.append(f"{r[0]} {r[1]} {r[2]} " + (f"{r[4]} {r[3]}" if k % 2 else f"{r[3]} {r[4]}"))
open(f"{w}/list.z", "w").write("\n".join(t) + "\n")
E
"$BIN" query --band imp.w100000000.t1.bld --out imp.list.ld --snps list.ids --precision 17 >/dev/null 2>&1
if out=$($CHK $W/imp imp.list.ld 100000000 0); then ok "SNP list (IDs): ${out#*: }"; else echo "$out"; bad "SNP list (IDs)"; fi
"$BIN" query --band imp.w100000000.t1.bld --out imp.listz.ld --snps list.z --precision 17 >/dev/null 2>&1
if out=$($CHK $W/imp imp.listz.ld 100000000 0) && [ "$(awk '$5!="G"' imp.listz.ld.snps.tsv | grep -vc '^#')" = 20 ]; then ok "SNP table with swapped alleles (signs flipped): ${out#*: }"; else echo "$out"; bad "SNP table with alleles"; fi

# filters and errors
gen mono 400 40 25 MONO=6 MISSFIX=9:3
"$BIN" build --bgen mono.bgen --out mono.b --window-bp 100000000 >mono.log 2>&1
[ "$(grep -vc '^#' mono.b.skipped.tsv)" = 2 ] && grep -q monomorphic mono.b.skipped.tsv && grep -q missing mono.b.skipped.tsv && ok "monomorphic and missing-data variants are skipped and listed" || bad "skip list: $(cat mono.b.skipped.tsv)"
"$BIN" info --band mono.b.bld | grep -q "Variants *: 38" && ok "info" || bad "info: $("$BIN" info --band mono.b.bld | tr '\n' ' ')"
head -c 3000 imp.w4000.t1.bld > trunc.bld
"$BIN" query --band trunc.bld --out t.ld --min-bp 0 --max-bp 2000000000 >/dev/null 2>err.txt && bad "truncated file accepted" || ok "truncated file: $(head -1 err.txt)"
"$BIN" query --band imp.w4000.t1.bld --out t.ld --snps /dev/null >/dev/null 2>err.txt && bad "empty list accepted" || ok "empty SNP list: $(head -1 err.txt)"

echo
echo "bandLD tests: $pass passed, $fail failed"
[ $fail -eq 0 ]
