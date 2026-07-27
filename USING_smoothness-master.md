# Building and using `localmethod` from Fortran source (any operating system)

`cp_source.R` calls compiled Fortran through `.Fortran(...)`, so before you can
`source()` it you must load a shared library built from the two Fortran files
(`localrisk.f90`, `localfit.f90`). That library is **machine-specific** — it is
tied to one operating system *and* one processor. A copy built elsewhere
(different OS, or Intel vs. Apple Silicon) will not load. The rule is simple:

> **Ship the `.f90` source, not the compiled binary. Rebuild on each machine.**

The build procedure below is the same on every OS; only the compiler you install
once differs.

---

## 1. Files you need

In one folder:

```
localrisk.f90     <- Fortran source
localfit.f90      <- Fortran source
cp_source.R       <- the R code that loads the library and calls it
```

Delete any pre-built binaries that came with the code before rebuilding, so a
stale one is never reused:

```
localmethod64.so   localmethod64.dll   *.o
```

(`.so` is the Unix/macOS name, `.dll` the Windows name — you only ever have one.)

---

## 2. One-time: install a Fortran compiler (the only OS-specific step)

A C compiler is **not** enough here — you need **gfortran** specifically.

- **Windows** — install **Rtools** for your R version
  (https://cran.r-project.org/bin/windows/Rtools/); it bundles gfortran. Accept
  the option to put it on the build PATH.
- **macOS** — install the **official gfortran** from
  https://mac.r-project.org/tools/ (choose the build matching your R version).
  It installs to `/opt/gfortran`, which is where R's configuration expects it.
  Xcode's Command Line Tools do **not** include Fortran.
- **Linux** — install via the package manager, e.g.
  Debian/Ubuntu: `sudo apt-get install gfortran r-base-dev`;
  Fedora/RHEL: `sudo dnf install gcc-gfortran R-devel`.

Check that R can find the compiler (works on every OS):

```
R CMD config FC        # should print a gfortran path that actually exists
```

If it prints a path that does **not** exist on your machine (common on macOS if
you installed gfortran via Homebrew instead of the official installer), point R
at the one you have by creating `~/.R/Makevars`:

```
FC  = /path/to/your/gfortran
F77 = /path/to/your/gfortran
FLIBS = -L/path/to/gfortran/lib -lgfortran -lquadmath -lm
```

(The lib path is `dirname $(gfortran -print-file-name=libgfortran.dylib)` on
macOS, or `...libgfortran.so` on Linux.)

---

## 3. Build the library (same command everywhere)

From the folder with the `.f90` files:

```
R CMD SHLIB -O3 localrisk.f90 localfit.f90 -o localmethod64.so
```

Let **R** drive the build — it invokes gfortran with the correct flags and
produces a library your R can load, with the right file type for your platform.
Notes that keep it portable:

- Keep the output name `localmethod64.so`. `cp_source.R` loads exactly that name
  (`dyn.load(file.path(SMOOTH_DIR, "localmethod64.so"))`), and R happily loads a
  `.so` on Windows too, so one name works on all systems and no R code changes.
- Do **not** add `--arch=...`; the default builds for the machine you are on,
  which is what you want.
- The code uses no OpenMP, so `-fopenmp` is unnecessary.
- On macOS, `-dynamiclib` is optional; `R CMD SHLIB` links a loadable library
  correctly either way.

This step recreates `localmethod64.so` (plus fresh `.o` files) for the current
machine.

---

## 4. Load and use it (identical everywhere)

In R, set `SMOOTH_DIR` to the folder holding the files, then:

```r
dyn.load(file.path(SMOOTH_DIR, "localmethod64.so"))
source(file.path(SMOOTH_DIR, "cp_source.R"))
```

---

## 5. Verify it loaded correctly

```r
dyn.load(file.path(SMOOTH_DIR, "localmethod64.so"))
is.loaded("localfit")        # TRUE
is.loaded("localrisk1")      # TRUE
```

All the routines `cp_source.R` relies on:
`localfit`, `localrisk1`, `localrisk2`, `iatra2`,
`localriskall`, `localriskall2`, `movingrisk`, `movingrisk2`, `movingrisk3`.

You can also confirm the architecture of the file you built:

```
file localmethod64.so        # macOS/Linux: shows arm64 / x86_64 / etc.
```

---

## 6. Troubleshooting

- **"incompatible architecture" or won't `dyn.load`** — the library was built for
  a different machine. Delete `localmethod64.so` and the `.o` files and redo
  step 3 on this machine.
- **`R CMD SHLIB`: "gfortran: No such file or directory"** — R points at a
  compiler you don't have; fix `R CMD config FC` via `~/.R/Makevars` (step 2).
- **`dyn.load` error "Library not loaded: libgfortran…"** — the runtime library
  can't be found; set `FLIBS` in `~/.R/Makevars` to the folder containing your
  `libgfortran` (step 2).
- **The Fortran source itself needs no edits** — it compiles cleanly under modern
  gfortran as-is.

---

### One-line summary

> Keep only the `.f90` source, install gfortran once for your OS, then run
> `R CMD SHLIB -O3 localrisk.f90 localfit.f90 -o localmethod64.so` and
> `dyn.load()` the result — R builds the correct library for whatever machine
> it runs on.
