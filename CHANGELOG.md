# Changelog

## [0.12.0](https://github.com/GeoGenetics/ngs-taxon/compare/v0.11.2...v0.12.0) (2026-08-24)


### Features

* implement dragen mapper ([#48](https://github.com/GeoGenetics/ngs-taxon/issues/48)) ([14e7948](https://github.com/GeoGenetics/ngs-taxon/commit/14e79484db3d8b37ac4eb4448c3aba5737a009cb))
* Implement onebam merge ([#53](https://github.com/GeoGenetics/ngs-taxon/issues/53)) ([51c268b](https://github.com/GeoGenetics/ngs-taxon/commit/51c268b72742a654bc8d658e4b8c715a80cfc709))
* Implement unicorn taxstats ([#55](https://github.com/GeoGenetics/ngs-taxon/issues/55)) ([5edd9ab](https://github.com/GeoGenetics/ngs-taxon/commit/5edd9ab63bdf1f37c014cca620d17582579c923f))
* remove bam_filter rules ([#54](https://github.com/GeoGenetics/ngs-taxon/issues/54)) ([29ce555](https://github.com/GeoGenetics/ngs-taxon/commit/29ce555db67923bf66eb5b7ea97600a7b64001c0))


### Bug Fixes

* dragen output paths ([#51](https://github.com/GeoGenetics/ngs-taxon/issues/51)) ([b8488b1](https://github.com/GeoGenetics/ngs-taxon/commit/b8488b1bb28513979a0cc2c52db9245d4ebb9b0b))


### Performance Improvements

* Bump Dragen ([#52](https://github.com/GeoGenetics/ngs-taxon/issues/52)) ([2b167e2](https://github.com/GeoGenetics/ngs-taxon/commit/2b167e28f1584d204b9fd2a84431e54be3124a44))
* bump metadmg ([#50](https://github.com/GeoGenetics/ngs-taxon/issues/50)) ([da1972b](https://github.com/GeoGenetics/ngs-taxon/commit/da1972b33f212afc73897eebe0aed47775f52428))

## [0.11.2](https://github.com/GeoGenetics/ngs-taxon/compare/v0.11.1...v0.11.2) (2026-03-06)


### Bug Fixes

* MultiQC benchmark extension ([#44](https://github.com/GeoGenetics/ngs-taxon/issues/44)) ([5650952](https://github.com/GeoGenetics/ngs-taxon/commit/565095287001cbc90feaded064049ddba50021f3))
* sort saturated reads to ensure passes unit testing ([#41](https://github.com/GeoGenetics/ngs-taxon/issues/41)) ([281c27e](https://github.com/GeoGenetics/ngs-taxon/commit/281c27e82127e5aa23a72f58f6a10a9a4e8ebfe8))


### Performance Improvements

* bump wrappers ([#46](https://github.com/GeoGenetics/ngs-taxon/issues/46)) ([82dd4a4](https://github.com/GeoGenetics/ngs-taxon/commit/82dd4a489ab22a01904ef0e32833b953eb40d4b1))
* update MultiQC version ([#45](https://github.com/GeoGenetics/ngs-taxon/issues/45)) ([2610adc](https://github.com/GeoGenetics/ngs-taxon/commit/2610adc1ea8c8d96e4a24533133aa5d4486c2e56))
* Update unicorn and env, and tweak resources and unit tests ([#42](https://github.com/GeoGenetics/ngs-taxon/issues/42)) ([b547900](https://github.com/GeoGenetics/ngs-taxon/commit/b547900c7fe07bd633021539074c440030d336f9))

## [0.11.1](https://github.com/GeoGenetics/ngs-taxon/compare/v0.11.0...v0.11.1) (2025-11-13)


### Bug Fixes

* Update pixi env and bump wrapper ver ([#39](https://github.com/GeoGenetics/ngs-taxon/issues/39)) ([613c07d](https://github.com/GeoGenetics/ngs-taxon/commit/613c07d223aa51f2c9278cc2fa56ee4c2aacfbdb))

## [0.11.0](https://github.com/GeoGenetics/ngs-taxon/compare/v0.10.1...v0.11.0) (2025-11-12)


### Features

* add bowtie2 env module to support CPU-specific binaries ([#33](https://github.com/GeoGenetics/ngs-taxon/issues/33)) ([c4b443f](https://github.com/GeoGenetics/ngs-taxon/commit/c4b443f36f81a16790fb7da5362c6873c986c568))


### Bug Fixes

* add missing samtools dep, and fix typo ([#36](https://github.com/GeoGenetics/ngs-taxon/issues/36)) ([08273b9](https://github.com/GeoGenetics/ngs-taxon/commit/08273b963258039e7961211aa87f0fa0ff99f0b5))


### Performance Improvements

* adjust resources, and bowtie2 exec options ([#35](https://github.com/GeoGenetics/ngs-taxon/issues/35)) ([ceb9daf](https://github.com/GeoGenetics/ngs-taxon/commit/ceb9daf106697a621cbce3a6e15f574288c24b5f))
* Reduce number of threads and bind to socket ([#38](https://github.com/GeoGenetics/ngs-taxon/issues/38)) ([1716446](https://github.com/GeoGenetics/ngs-taxon/commit/1716446d1c6d9ccaa7bdc2053407b557ee778728))

## [0.10.1](https://github.com/GeoGenetics/ngs-taxon/compare/v0.10.0...v0.10.1) (2025-07-25)


### Bug Fixes

* revert metadmg ([#31](https://github.com/GeoGenetics/ngs-taxon/issues/31)) ([e6059f5](https://github.com/GeoGenetics/ngs-taxon/commit/e6059f5a5812f7578fd8de8a156ef5213c58454b))

## [0.10.0](https://github.com/GeoGenetics/ngs-taxon/compare/v0.9.6...v0.10.0) (2025-07-24)


### Features

* Bump metadmg version ([#30](https://github.com/GeoGenetics/ngs-taxon/issues/30)) ([0767200](https://github.com/GeoGenetics/ngs-taxon/commit/0767200c60949e808cf4ac1197e4f024de1d8ccd))
* replace compressbam with unicorn ([#27](https://github.com/GeoGenetics/ngs-taxon/issues/27)) ([7e63ed1](https://github.com/GeoGenetics/ngs-taxon/commit/7e63ed1d98756e769438dcc633ebfb3a64a65c48))

## [0.9.6](https://github.com/GeoGenetics/ngs-taxon/compare/v0.9.5...v0.9.6) (2025-07-19)


### Bug Fixes

* tweak metadmg_dfit resources ([#25](https://github.com/GeoGenetics/ngs-taxon/issues/25)) ([c80c65e](https://github.com/GeoGenetics/ngs-taxon/commit/c80c65e9e3bb603ad2ec0819ebf9383252e0d2a1))

## [0.9.5](https://github.com/GeoGenetics/ngs-taxon/compare/v0.9.4...v0.9.5) (2025-07-18)


### Bug Fixes

* resource typo ([#23](https://github.com/GeoGenetics/ngs-taxon/issues/23)) ([4e9e45f](https://github.com/GeoGenetics/ngs-taxon/commit/4e9e45f92894005207a1218b3b24cdb3aec902d3))

## [0.9.4](https://github.com/GeoGenetics/ngs-taxon/compare/v0.9.3...v0.9.4) (2025-07-18)


### Bug Fixes

* default test units ([#18](https://github.com/GeoGenetics/ngs-taxon/issues/18)) ([08f28d2](https://github.com/GeoGenetics/ngs-taxon/commit/08f28d240ec22c4c385e6020dbc233a62bdf76d5))
* minimum runtime ([#21](https://github.com/GeoGenetics/ngs-taxon/issues/21)) ([a304812](https://github.com/GeoGenetics/ngs-taxon/commit/a3048126da6bdc49e39e3fa863a5a6a0141cec21))
* saturated reads output ([#20](https://github.com/GeoGenetics/ngs-taxon/issues/20)) ([1339b7c](https://github.com/GeoGenetics/ngs-taxon/commit/1339b7ca6d0a5d0a8ebf7ab62ef844297ae16aab))

## [0.9.3](https://github.com/GeoGenetics/ngs-taxon/compare/v0.9.2...v0.9.3) (2025-06-30)


### Bug Fixes

* bam_filer query sorting ([#15](https://github.com/GeoGenetics/ngs-taxon/issues/15)) ([a922244](https://github.com/GeoGenetics/ngs-taxon/commit/a92224443bebe7bbd63fb52110d95a3f7dc2bd55))
* use metadmg conda  ([#17](https://github.com/GeoGenetics/ngs-taxon/issues/17)) ([ee748f4](https://github.com/GeoGenetics/ngs-taxon/commit/ee748f49ff10bc7e99d5fc9b21bf5963678aad1b))

## [0.9.2](https://github.com/GeoGenetics/ngs-taxon/compare/v0.9.1...v0.9.2) (2025-06-19)


### Bug Fixes

* increase bam_filter resources and comment param ([#13](https://github.com/GeoGenetics/ngs-taxon/issues/13)) ([32a00fe](https://github.com/GeoGenetics/ngs-taxon/commit/32a00fe8ca72401a52e13c6f7d6a6ce6665764e5))

## [0.9.1](https://github.com/GeoGenetics/ngs-taxon/compare/v0.9.0...v0.9.1) (2025-06-19)


### Bug Fixes

* validation URL and schema ([#11](https://github.com/GeoGenetics/ngs-taxon/issues/11)) ([f9a81d5](https://github.com/GeoGenetics/ngs-taxon/commit/f9a81d547def9fa1e8dc7be2a957e0d7a1f0c4ca))

## [0.9.0](https://github.com/GeoGenetics/ngs-taxon/compare/v0.8.1...v0.9.0) (2025-06-17)


### Features

* Add filtering of saturated reads ([#7](https://github.com/GeoGenetics/ngs-taxon/issues/7)) ([2294dc3](https://github.com/GeoGenetics/ngs-taxon/commit/2294dc3dea3804714824db91e25d7844b539886d))


### Bug Fixes

* Code refactor, tweak resources, bump wrappers, and update pixi env. ([#8](https://github.com/GeoGenetics/ngs-taxon/issues/8)) ([bdff102](https://github.com/GeoGenetics/ngs-taxon/commit/bdff1026aa14f747803e134177fa2ba6cbb5f71a))
* YAML schema ([#9](https://github.com/GeoGenetics/ngs-taxon/issues/9)) ([1e5adc3](https://github.com/GeoGenetics/ngs-taxon/commit/1e5adc31e8a286d802f87bc27ebfea02c554fa8f))


### Performance Improvements

* remove slurm nice for mapping ([#6](https://github.com/GeoGenetics/ngs-taxon/issues/6)) ([f205f7d](https://github.com/GeoGenetics/ngs-taxon/commit/f205f7d5899f489cf7254a68cc24fec064ad28f7))
* shard sort query ([#5](https://github.com/GeoGenetics/ngs-taxon/issues/5)) ([7789abe](https://github.com/GeoGenetics/ngs-taxon/commit/7789abef488e5f6463bed2ae6c1be41fc2ed2302))
* Tweak resources ([#2](https://github.com/GeoGenetics/ngs-taxon/issues/2)) ([25c0b8d](https://github.com/GeoGenetics/ngs-taxon/commit/25c0b8d6868c72e353db8efa5b34e93fd9fe0096))
