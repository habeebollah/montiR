# montiR
montiR (Model biOmass diNamik singkaT dI R) merupakan kumpulan model biomas dinamik dengan pendekatan time series fitting. Tool pemodelan yang tersedia disini antara lain:

1. Data plotting

Langkah paling penting sebelum melakukan analisis data adalah memeriksa apakah data yang akan digunakan memenuhi persyaratan dan asumsi yang dibutuhkan untuk analisis biomass dynamic model, termasuk memilih jenis langkah apa yang harus dilakukan ketika data yang dibutuhkan tidak memenuhi asumsi

2. Surplus produksi dengan asumsi equilibrium

Tool ini menghitung jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) dengan asumsi equilibrium untuk model Schaefer dan Fox. Pendekatan ini ditampilkan disini hanya untuk tujuan edukasi sebagai contoh model yang akan memberikan estimasi MSY dan Emsy yang lebih tinggi, sehingga sangat tidak disarankan untuk dijadikan sebagai panduan dalam pengambilan kebijakan perikanan. Overestimasi reference point pada kondisi ketika status perikanan sedang dalam kondisi overexploited akan memberikan ilusi bahwa stok ikan masih banyak, sehingga dapat merugikan pelaku perikanan karena jumlah tangkapan yang rendah dan merugikan stok ikan karena semakin tingginya pemanfaatan. Banyak dari kita yang masih menggunakan metode ini, meskipun sudah tidak disarankan untuk digunakan sejak 1980an.

3. Surplus produksi dengan asumsi non equilibrium menggunakan metode time series fitting

Tool ini melakukan estimasi parameter K, B0, r, q dan menentukan jumlah tangkapan ikan lestari (MSY), biomassa ikan lestari (Bmsy), serta upaya penangkapan ikan lestari (Emsy) menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox. Tool ini sudah disesuaikan untuk kebutuhan data yang terbatas (dapat mengakomodasi ketiadaaan input data upaya penangkapan) serta sudah memperhitungkan kesalahan dalam pengambilan data (observation error) dan kesalahan dalam model (model error). Metode time series fitting disebut sebagai metode yang lebih baik dibandingkan dengan dua metode lain (metode equilibrium dan regresi) yang digunakan untuk melakukan estimasi parameter dalam model surplus produksi. Sebagian kecil dari kita sudah menggunakan metode ini, tetapi masih kurang tepat dalam melakukan analisisnya sehingga berakibat pada kurang tepatnya perhitungan MSY, Bmsy dan Emsy

4. Menghitung rentang data reference point

Tool ini melakukan estimasi atas reference point untuk jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) serta menghitung standard error menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox. Pendekatan ini akan memberikan rentang kemungkinan error atas pendugaan MSY dan Emsy yang dilakukan.

5. Menyediakan data prior untuk parameter pertumbuhan r

Data prior setiap spesies untuk parameter pertumbuhan r diambil dari database fishbase dan sealifebase untuk mendukung pendugaan stok di tingkat spesies

6. Menghitung reference point untuk pengelolaan di tingkat spesies

Metode ini digunakan untuk menentukan MSY, Bmsy dan Emsy menggunakan pendekatan bayesian dan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox sehingga menghasilkan pendugaan stok yang lebih akurat di tingkat spesies

7. Membuat proyeksi atas kebijakan reference point berdasar tingkat pemanfaatan perikanan

Tool ini akan membuat grafik proyeksi biomass per biomass at msy (B/Bmsy) dan fishing per fishing at msy (F/Fmsy) sebagai panduan untuk melihat kebijakan yang dibuat saat ini serta memperkirakan limit dan target reference point

## Installation

You can install the development version of montiR like so:

``` r
devtools::install_github("habeebollah/montiR")
```

## Example

This is a basic example which shows you how to solve a common problem:

```{r example}
library(montiR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

```{r cars}
summary(cars)
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>.
