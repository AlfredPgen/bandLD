# bandLD

Chromosome-wide LD in one compact file. `bandLD build` stores the LD of every pair of variants that
are closer than a window (default 1 Mb) along a whole chromosome. `bandLD query` returns the full
LD matrix of any region up to that width, or of any list of variants.

Version 0.1.0, first working version. Linux x86-64. Input: BGEN v1.2/1.3 with a bgenix `.bgi`
index (dosages or hard calls, biallelic, diploid).

## Install

Download `bandLD_linux_x86_64` from the [releases page](https://github.com/AlfredPgen/bandLD/releases)
and run `chmod +x bandLD_linux_x86_64`. It is a static binary with no dependencies and picks
AVX2 or AVX-512 kernels at run time. To compile instead (Debian/Ubuntu):

```
sudo apt install g++ make libdeflate-dev libzstd-dev libsqlite3-dev
make            # or: make static
```

## Build a band file

```
bandLD build --bgen chr2.bgen --out chr2 --threads 16
```

Writes `chr2.bld` (the band file), `chr2.snps.tsv` (one line per stored variant: chromosome,
position, rsid, allele1, allele2, MAF, row index) and `chr2.skipped.tsv` (variants left out, with
the reason). The LD sign always refers to `allele2`. To index the table:

```
bgzip chr2.snps.tsv && tabix -s 1 -b 2 -e 2 chr2.snps.tsv.gz
```

The run prints its settings, the expected memory, and a progress line every 30 seconds (variants
done, pairs written, file size and bytes per pair, pairs per second, time elapsed and left).

Options: `--window-bp N` (default 1000000), `--bgi FILE`, `--chr NAME`, `--sample FILE` and
`--incl FILE` (sample subset), `--extract FILE` (variant IDs), `--min-maf X`, `--min-bp N --max-bp N`
(part of a chromosome), `--round-bits B` and `--exact` (see Precision).

## Query

```
bandLD query --band chr2.bld --out region.ld.gz --chr 2 --min-bp 12345 --max-bp 1012345
bandLD query --band chr2.bld --out mylist.ld --snps mylist.txt
bandLD info  --band chr2.bld
```

The matrix is written as text (gzipped if the name ends in `.gz`), in the same layout as fastLD,
with `OUT.snps.tsv` listing the variants in matrix order. `--snps` takes variant IDs, one per
line, or a table with a header line containing `rsid position allele1 allele2`; with the table,
alleles are checked and the sign follows the table's `allele2`. Pairs further apart than the window
are not stored and print as `nan` (a warning gives the count).

## Precision

* **16-bit (default).** r is rounded to a step of 2^-15 (3e-5), the precision of 16-bit storage.
  The largest error, 1.5e-5, is about 200 times smaller than the sampling error of r at 100,000
  samples. `--round-bits B` chooses another step of 2^-B.
* **Exact (`--exact`).** For each pair the file holds the exact integer cross product of the two
  variants' dosages; with the per-variant sums in the file this gives back r as the same double that
  fastLD computes. A region queried from an exact band file is byte-identical to fastLD's output on
  the same BGEN.

## Size and speed (measured)

Simulated 1 Mb region, 9,345 variants, 100,000 samples, 43.7 million pairs, laptop CPU
(i7-11370H). Bytes per pair:

| Precision | bandLD, dosages | bandLD, hard calls | Zarr method, dosages |
|---|---|---|---|
| exact (`--exact`) | 2.83 | 0.80 | not available |
| step 2^-20 (`--round-bits 20`) | 1.95 | 1.88 | |
| 16-bit (default) | 1.33 | 1.25 | 1.44 |
| step 2^-11 (`--round-bits 11`) | 0.83 | 0.76 | |
| 8-bit (`--round-bits 7`) | 0.35 | 0.34 | 0.43 |
| step 2^-3 (`--round-bits 3`) | 0.07 | 0.08 | |

"Zarr method" = the same values rounded to int16 or int8 and compressed with Blosc (zstd level 7,
byte shuffle), which is what magenpy writes. The simulated dosage noise is a guess; run a pilot on
real data (`--min-bp/--max-bp` over 2 Mb) to get the true figure before a full build.

* Build: 95 s on one thread, 35 s on four (about 510,000 pairs per second per thread at this sample
  size); memory 2.1 GB. Output is identical for any thread count and CPU kernel.
* Query: a 200 kb region (1,886 variants) in 1.3 s, a list of 100 variants in 0.9 s, the full 1 Mb
  region (9,345 variants) in 11 s.

## How it works

Variants are read in position order. A sliding store keeps the dosages of every variant still within
the window of the current position (2 bytes per sample per variant), and each new variant's exact
cross products with all of them are computed with fastLD's integer kernels (AVX-512 VNNI, AVX2 or
scalar). Row j of the file holds the values for all earlier variants in its window. Values are
entropy coded with an adaptive range coder whose statistics depend on the expected noise level of
each pair and on the distance between the two variants; rows are grouped in chunks of 64 so a query
decodes only the chunks it needs. Format details: `src/bandfile.hpp`.

## Limits of this version

* Variants with missing genotypes or with probabilities of more than 10 bits are skipped (listed in
  `.skipped.tsv`), as are monomorphic variants.
* Only the dense kernel is used; low-frequency variants are not yet stored sparsely.
* BGEN input only. No R reader yet (the query output is the text matrix).

## Tests

`make test` (python3, numpy, scipy): every stored value is compared with an exact reference for
several window sizes, thread counts and kernels; rounded mode; region and list queries; allele
swaps; skipped variants; corrupt files. `make asan` builds with AddressSanitizer and UBSan.

BGEN decoding and the integer kernels are taken from fastLD (https://github.com/AlfredPgen/fastLD).

## Licence

MIT, see `LICENSE`. Libraries linked into the release binary are listed in `THIRD_PARTY.md`.
