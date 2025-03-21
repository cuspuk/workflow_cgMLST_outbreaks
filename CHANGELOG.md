# Changelog

## [1.6.3](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.6.2...v1.6.3) (2025-03-21)


### Bug Fixes

* added neisseria into config ([cc8d9a8](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/cc8d9a8dae8e2a7c18b4dab730f3d56893aa82d2))
* added neisseria into config ([13f05b6](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/13f05b664b6d13bb31ee27455ed1e35237bf797a))

## [1.6.2](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.6.1...v1.6.2) (2024-11-15)


### Bug Fixes

* case when samples less than 2 ([e164936](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/e164936f7e8ee8cb3cb0ec2c30aa5168db4a05ab))

## [1.6.1](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.6.0...v1.6.1) (2024-11-12)


### Bug Fixes

* moved samples ini into metadata dir ([a286fba](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/a286fba5d6557dd6eb98daaca8cd9dc680015da0))
* moved samples under cgmlst dir ([362bdea](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/362bdea213b4ac42298ef7512b922bff414631e7))


### Performance Improvements

* marked dump sampling as localrule ([951cb3f](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/951cb3fa2af4315acf97759597585305200a49da))

## [1.6.0](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.5.0...v1.6.0) (2024-11-12)


### Features

* added logging of samples used per each taxa ([8406c57](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/8406c574a2ec943637dd9552457780a3d9eb389b))


### Bug Fixes

* newick env ([1da4877](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/1da4877c23708c0259f73b99a801213e1499852c))

## [1.5.0](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.4.2...v1.5.0) (2024-11-11)


### Features

* added newick tree ([f85ce83](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/f85ce83cf435ed5a0b98cb977c2b52365ca59621))

## [1.4.2](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.4.1...v1.4.2) (2024-07-27)


### Bug Fixes

* bumped chewbacca version ([9ea07ae](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/9ea07ae56ad46c75be228d44e60483089f5deb69))

## [1.4.1](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.4.0...v1.4.1) (2024-07-25)


### Bug Fixes

* prevent redownloading and removing custom paths ([f5c100c](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/f5c100c8d9cdb8fd31ad4a15198ded8efe4de576))

## [1.4.0](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.3.1...v1.4.0) (2024-07-13)


### Features

* updated config and coupled trn files ([2c43bae](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/2c43baec05e5a64d74ca7f5b5c2e27bd10c9ba4a))

## [1.3.1](https://github.com/cuspuk/workflow_cgMLST_outbreaks/compare/v1.3.0...v1.3.1) (2024-07-08)


### Bug Fixes

* added checks for urls if schemas not already present. Added check for prefix ([074bd29](https://github.com/cuspuk/workflow_cgMLST_outbreaks/commit/074bd29e602fc89ffdfdd1af143bd37bd9c260f7))

## [1.3.0](https://github.com/xsitarcik/cgMLST_outbreaks/compare/v1.2.1...v1.3.0) (2024-03-27)


### Features

* refactored config to accept full paths instead of relative paths for schemas ([6ad4117](https://github.com/xsitarcik/cgMLST_outbreaks/commit/6ad41174c844cf402a6b360cfb94a7e7f7f1ab5c))


### Bug Fixes

* request outputs only for taxa labels that have samples ([bc7de9e](https://github.com/xsitarcik/cgMLST_outbreaks/commit/bc7de9edf106a663b435cffc9e0ad48a958eb268))

## [1.2.1](https://github.com/xsitarcik/cgMLST_outbreaks/compare/v1.2.0...v1.2.1) (2024-03-07)


### Bug Fixes

* fix chewbacca handling of bad filenames ([c28afd1](https://github.com/xsitarcik/cgMLST_outbreaks/commit/c28afd1de7dd87c839ddccf52216dc184fb4e800))
* pepfile script works with analyses from assembly ([6be4c52](https://github.com/xsitarcik/cgMLST_outbreaks/commit/6be4c523f3ade009a1328b634fe20cf46ed8b0be))

## [1.2.0](https://github.com/xsitarcik/cgMLST_outbreaks/compare/v1.1.0...v1.2.0) (2024-03-07)


### Features

* added cgmlst_dists ([a985162](https://github.com/xsitarcik/cgMLST_outbreaks/commit/a98516225a06c16f32a145dc42a49f7045371c3f))
* reworked to allow other schemas than ridom ([69e99c2](https://github.com/xsitarcik/cgMLST_outbreaks/commit/69e99c2c109d1976953429ff55e72e7d628b57dd))
* samples are processed per taxa list ([ff166e9](https://github.com/xsitarcik/cgMLST_outbreaks/commit/ff166e99f266ba8058a74a3185a88adeb80fbc0c))

## [1.1.0](https://github.com/xsitarcik/cgMLST_outbreaks/compare/v1.0.0...v1.1.0) (2024-03-06)


### Features

* added script for building pepfile ([5638730](https://github.com/xsitarcik/cgMLST_outbreaks/commit/56387301f787b0fc6c5489ff930f60778c8b4e69))

## 1.0.0 (2024-03-05)


### Features

* added init version ([7c333b9](https://github.com/xsitarcik/cgmlst_outbreaks/commit/7c333b9c6064fd00aa9b822cc101ea28850e0375))


### Bug Fixes

* chewbacca output removed before calling to ensure clean results ([b1a8d15](https://github.com/xsitarcik/cgmlst_outbreaks/commit/b1a8d155dbd624538df1fdac1a96c1657ae73fc2))
