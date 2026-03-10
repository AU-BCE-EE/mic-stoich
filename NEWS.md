# micstoich 0.1.1

* Fixed bug in internal `massconv()` function where `calcCOD()` was called without `per = 'mol'`, causing incorrect unit conversion for compounds with COD' > 0. This is a non-exported function.
