# DiscovarDenovo
This repo contains a copy of the DiscovarDenovo short-read assembler developed by Broad Institute.
The repo is a mirror of [discovarexp-51885](ftp://ftp.broadinstitute.org/pub/crd/DiscovarDeNovo/latest_source_code/discovarexp/discovarexp-51885.tar.gz)
 
## Note for building process
Builds with modern compilers (tested with GCC 16 under `-std=gnu++14`):

    autoreconf -fi && ./configure && make -j

jemalloc is auto-detected when present (`--with-jemalloc=DIR` to force a
location); the legacy GCC 4.8 requirement and its dead compiler flags are
gone. Developer notes: `docs/GraphConstruction.md` (graph construction +
memory optimization, items M1/M2 now implemented), `docs/Modernization.md`
(roadmap; P0 toolchain/memory work landed).

DISCOVAR de novo (Experimental)
===============================

DISCOVAR de novo is a short read genome assembly developed by the Broad Institute. It can generate high quality assemblies from the latest 250 base Illumina paired-end PCR-free fragment reads. This is an early release and is not yet completely stable or fully functional, hence the Experimental tag.

DISCOVAR is undergoing continuous development. There is no 'official' release, instead please download the latest source code from our nightly builds. Check back often for updates.


The LICENSE file contains copyright and licensing information.

The INSTALL file contains information on how to configure, build, and install DISCOVAR de novo.


The docs directory contains documentation on DISCOVAR de novo and its associated tools:

    DiscovarExp      - DISCOVAR de novo manual
    NhoodInfo        - NhoodInfo manual
    CrossOut         - CrossOut manual
    GraphConstruction.md - developer notes: unitig graph construction (K=60 -> K=200)
                       and memory optimization of the a.200 stage
    Modernization.md  - developer roadmap: shallow long-read hybrid assembly
                       (ONT/PacBio) and phasing on top of the Discovar
                       short-read core

Developer documentation on how short-read error correction is performed
inside the codebase is in [docs/ErrorCorrection.md](docs/ErrorCorrection.md).


For more information, visit the DISCOVAR blog:
http://www.broadinstitute.org/software/discovar/blog/

The above link includes directions on using the tool, along with other documentation.


-Computational R&D, The Broad Institute, Cambridge, MA
http://www.broadinstitute.org/science/programs/genome-biology/computational-rd/computational-research-and-development

