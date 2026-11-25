# Beagle JAR

This folder holds `beagle.27Feb25.75f.jar`, Beagle 5.5, the program that `cgiarGenomics::impute_beagle()` runs. The function finds it through `system.file("extdata", package = "cgiarGenomics")`. To use another release, place its JAR here or pass its path as `beagle_path`.

| | |
|---|---|
| Program | Beagle 5.5, version 27Feb25.75f |
| Author | Brian L. Browning |
| Source and download | https://faculty.washington.edu/browning/beagle/beagle.html |
| SHA-256 of the file in this folder | `7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70` |

## Licence

Beagle is free software under the GNU General Public License, version 3 or any later version (https://www.gnu.org/licenses/gpl-3.0.html). It is a separate program. cgiarGenomics starts it through Java and does not link to it. The MIT licence of cgiarGenomics does not apply to this file.

Anyone who redistributes this folder passes on the JAR under the GPL, with the licence text and access to the Beagle source code, which the page above provides.

Cite Beagle when you publish results imputed with it. The page above lists the articles.
