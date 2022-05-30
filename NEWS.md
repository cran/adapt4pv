## version version 0.2-2

---

- Add reference to the DESCRIPTION file.

## version version 0.2-1

---

- Bug fixed in the internal function hdps_pv: the denominator of the PC1 quantity is Tot_E and not Tot.

## version version 0.2-0

---

- Modification on all functions that involves cross-validation with the adaptive lasso  (cf. Nadim, B., Lola, E., and Vivian, V. (2020). On the use of cross-validation for the calibration of the tuning parameter in the adaptive lasso. arXiv preprint arXiv:2005.10119.)
- The adapt_ridge function has been removed, the adapt_cv and adapt_univ (with criterion = cv) functions have been modified.
- In this version, two differents implementation of the cross-validation for th adaptive lasso are possible through the new parameter type_cv. When type_cv = naive, it corresponds to the first version of this package: adaptive weights obtained on the full data are used for the cross-validation. When type_cv = proper, adaptive weights are calculated for each training sets.


