## Test environments

(updated for v1.1.1)

-   Local: OS X install, R version 4.0.3
-   AppVeyor: <https://ci.appveyor.com/project/stineb/rpmodel>: Successful build of v1.1.1

## R CMD check results

### devtools::check() locally:

0 errors ✓ \| 0 warnings ✓ \| 0 notes ✓

### devtools::check_rhub()

'OK' for all three platforms.

0 errors ✓ \| 0 warnings ✓ \| 0 notes ✓

### devtools::check_win_devel()

Status: 1 NOTE

Found the following (possibly) invalid file URI: URI: ./articles/usage.html
Unclear why this is raised. It works in the compiled html (package website built with pkgdown::build_site).

## Downstream dependencies

-   rlang, version 0.4.0 successfully installed.
