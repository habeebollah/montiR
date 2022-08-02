# montiR
montiR (Model biOmass diNamik singkaT dI R) merupakan kumpulan model biomas dinamik dengan pendekatan time series fitting. Tool pemodelan yang tersedia disini antara lain:

(1) Surplus produksi dengan asumsi non-equilibrium menggunakan data fitting

(a) Data plotting

Langkah paling penting sebelum melakukan analisis data adalah memeriksa apakah data yang akan digunakan memenuhi persyaratan dan asumsi yang dibutuhkan untuk analisis biomass dynamic model, termasuk memilih jenis langkah apa yang harus dilakukan ketika data yang dibutuhkan tidak memenuhi asumsi

(b) Surplus produksi dengan asumsi non-equilibrium menggunakan time series fitting

Tool ini melakukan estimasi parameter K, B0, r, q dan menentukan jumlah tangkapan ikan lestari (MSY), biomassa ikan lestari (Bmsy), serta upaya penangkapan ikan lestari (Emsy) menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer. Tool ini sudah disesuaikan untuk kebutuhan data yang terbatas (dapat mengakomodasi hilangnya input data upaya penangkapan) serta sudah memperhitungkan kesalahan dalam pengambilan data (observation error) sehingga meningkatkan akurasi estimasi stok. Metode time series fitting disebut sebagai metode yang lebih baik dibandingkan dengan dua metode lain (metode equilibrium dan multiple regression) yang digunakan untuk melakukan estimasi parameter dalam model surplus produksi. Sebagian kecil dari kita sudah menggunakan metode ini, tetapi pelaksanaan analisisnya masih perlu disempurnakan agar tidak berakibat pada kurang tepatnya perhitungan MSY, Bmsy dan Emsy

(c) Menghitung reference point untuk pengelolaan

Tool ini menghasilkan jumlah stok ikan yang lestari (Bmsy), jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) untuk model Schaefer.

(d) Menghitung standard error dari reference point

Tool ini mengestimasi jumlah stok ikan yang lestari (Bmsy), jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) serta menghitung standard error menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer.

(2) Meningkatkan akurasi pendugaan stok dengan surplus production di tingkat spesies menggunakan bayesian

Data prior setiap spesies untuk parameter pertumbuhan r diambil dari database fishbase dan sealifebase untuk mendukung pendugaan stok di tingkat spesies. Langkah lanjutan untuk melakukan estimasi parameter surplus production menggunakan pendekatan bayesian dan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox sehingga menghasilkan pendugaan stok yang lebih akurat di tingkat spesies juga diberikan.

(3) Membuat proyeksi atas kebijakan reference point berdasar tingkat pemanfaatan perikanan

Tool ini akan membuat grafik proyeksi biomass per biomass at msy (B/Bmsy) dan fishing per fishing at msy (F/Fmsy) sebagai panduan untuk melihat kebijakan yang akan dibuat berdasarkan Bmsy, MSY dan Emsy sebagai reference point. Proyeksi dibuat dengan pendekatan deterministic secara default, dan terdapat opsi untuk tujuan stochastic.

## Installation

Anda dapat mengunduh dan menginstal package ini dengan mengetikkan kode ini pada konsole di Rstudio:

``` r
devtools::install_github("habeebollah/montiR")
```

## Panduan penggunaan

Panduan penggunaan package ini dapat dilihat pada link berikut:

``` r
https://habeebollah.github.io/
```
