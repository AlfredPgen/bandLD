# Third-party software

bandLD's sources are MIT-licensed (see `LICENSE`). The BGEN decoder and the integer kernels come
from fastLD (same author, MIT). The release binary is statically linked and contains:

| Library | Used for | Licence |
|---|---|---|
| libdeflate | zlib decompression of BGEN, gzip output | MIT |
| Zstandard | zstd decompression of BGEN, SNP table | BSD-3-Clause |
| SQLite | reading bgenix `.bgi` indexes | public domain |
| GCC runtime (libstdc++, libgomp, libgcc) | runtime | GPL-3.0 with the GCC Runtime Library Exception |
| GNU C Library | runtime | LGPL-2.1-or-later |

The complete source is in this repository, so the binary can be rebuilt and relinked with
`make static`.
