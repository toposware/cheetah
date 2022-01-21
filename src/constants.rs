// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::{BasePointTable, LookupTable, ProjectivePoint};
use crate::{Fp, Fp2, Fp6};

/// A hardcoded `BasePointTable` for the generator of the
/// Cheetah curve, to allow for efficient scalar multiplication.
#[doc(hidden)]
pub const BASEPOINT_TABLE: BasePointTable = BasePointTable([
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xf6798582c92ece1),
                    c1: Fp(0x2b7c30a4c7d886c0),
                },
                c1: Fp2 {
                    c0: Fp(0x1269cdae98dc2fd0),
                    c1: Fp(0x11b78ef6c71c6132),
                },
                c2: Fp2 {
                    c0: Fp(0x3ac2244dfc47537),
                    c1: Fp(0x36dfeea4b9051daf),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x334807e450d55e2f),
                    c1: Fp(0x200a54d42b84bd17),
                },
                c1: Fp2 {
                    c0: Fp(0x271af7bb20ab32e1),
                    c1: Fp(0x3df7b90927efc7ec),
                },
                c2: Fp2 {
                    c0: Fp(0xab8bbf4a53af6a0),
                    c1: Fp(0xe13dca26b2ac6ab),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ba2d52806f212a),
                    c1: Fp(0x5e9353a4e8225c8),
                },
                c1: Fp2 {
                    c0: Fp(0x13e92423fef3bc2d),
                    c1: Fp(0x241081e7ae1db310),
                },
                c2: Fp2 {
                    c0: Fp(0x29f0073c3351026b),
                    c1: Fp(0x11233fe9eb7285c0),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a19dfba18e15ed5),
                    c1: Fp(0x3691eb6949fca20b),
                },
                c1: Fp2 {
                    c0: Fp(0x3ea42cb9ad7430ab),
                    c1: Fp(0x1b840f91119a2eb3),
                },
                c2: Fp2 {
                    c0: Fp(0x1b94f8ccdafc47ba),
                    c1: Fp(0x19e92e12c3a9cfa),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x240fa5ba0a5898ed),
                    c1: Fp(0x40c47d7ba253bf0e),
                },
                c1: Fp2 {
                    c0: Fp(0x33ba4aded4efccc),
                    c1: Fp(0x1039b991fc425c2a),
                },
                c2: Fp2 {
                    c0: Fp(0x3a3377cb089e18e1),
                    c1: Fp(0x33fe6dd67a1199d1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3d18545ae8d7f1ca),
                    c1: Fp(0x35d016be897d1ac),
                },
                c1: Fp2 {
                    c0: Fp(0x11cffedc95fab7c),
                    c1: Fp(0x380a47d8d19fe89d),
                },
                c2: Fp2 {
                    c0: Fp(0xf2626400fe07952),
                    c1: Fp(0x29bdee7543df2dc1),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x38421f520be9eaa5),
                    c1: Fp(0x42de68de893c711),
                },
                c1: Fp2 {
                    c0: Fp(0x3e8e067736bd2eb8),
                    c1: Fp(0x30460eddb629aa0b),
                },
                c2: Fp2 {
                    c0: Fp(0x12d35376f4dda11c),
                    c1: Fp(0x3900acab8b415e26),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1d8b65a1650963c8),
                    c1: Fp(0x1a9e0c9571c0cd5e),
                },
                c1: Fp2 {
                    c0: Fp(0xffb279b803f9a1e),
                    c1: Fp(0x3f4362974ac3f3e4),
                },
                c2: Fp2 {
                    c0: Fp(0x3b1c67d6d589d25e),
                    c1: Fp(0xd16b2f32f47b2c5),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x20963a9466642969),
                    c1: Fp(0x6d053c1cc2028ae),
                },
                c1: Fp2 {
                    c0: Fp(0x20113cd3b01b2e42),
                    c1: Fp(0x3095357e8999801),
                },
                c2: Fp2 {
                    c0: Fp(0x16a416abf236f9e7),
                    c1: Fp(0x2f3a322b2d83abbd),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ef655cb5df8e10e),
                    c1: Fp(0x288c475d23a8a254),
                },
                c1: Fp2 {
                    c0: Fp(0x1dff0ce61a7ae2bc),
                    c1: Fp(0x1129ae857ee07f8a),
                },
                c2: Fp2 {
                    c0: Fp(0x870bd430eca1e97),
                    c1: Fp(0x1cd4733a0462ae02),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x21f56b2c2c2b9f74),
                    c1: Fp(0x26c1f4f022e38324),
                },
                c1: Fp2 {
                    c0: Fp(0x2fff141ac2c11242),
                    c1: Fp(0x18f66618417b4ff4),
                },
                c2: Fp2 {
                    c0: Fp(0x2aecbec97c6022fa),
                    c1: Fp(0x128adabb230d9814),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x4559c3329006d77),
                    c1: Fp(0xd50ac681186a14d),
                },
                c1: Fp2 {
                    c0: Fp(0x1dc79c586ab3a07),
                    c1: Fp(0x15c428abe161368f),
                },
                c2: Fp2 {
                    c0: Fp(0x4c6f6b905cf7c50),
                    c1: Fp(0x3fe15445d1cea5b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1779944f15ceab03),
                    c1: Fp(0x2e5b13cb48482f81),
                },
                c1: Fp2 {
                    c0: Fp(0x4073839e4414218d),
                    c1: Fp(0x1b99bde9130cd32a),
                },
                c2: Fp2 {
                    c0: Fp(0xa6a4605df9774e4),
                    c1: Fp(0x3d98693f8ca0016),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x22d756e15bca4b42),
                    c1: Fp(0x29703b4435c7ea5),
                },
                c1: Fp2 {
                    c0: Fp(0x399d43e5b876ab4),
                    c1: Fp(0x172962906cb0b347),
                },
                c2: Fp2 {
                    c0: Fp(0x21bcc1674db9f7ba),
                    c1: Fp(0x272437362226a1b4),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c6db9906b674248),
                    c1: Fp(0x93348e2de77d662),
                },
                c1: Fp2 {
                    c0: Fp(0x2b098441acf17bd1),
                    c1: Fp(0x2836d01fa21c3267),
                },
                c2: Fp2 {
                    c0: Fp(0x36a95afba527cbda),
                    c1: Fp(0x3170d4d25d1fc074),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3dad8cd910195fa4),
                    c1: Fp(0x18fda15837f0c55d),
                },
                c1: Fp2 {
                    c0: Fp(0x24a92bc9469374b8),
                    c1: Fp(0x1fc214214fa71c00),
                },
                c2: Fp2 {
                    c0: Fp(0x44b0f35a91ce992),
                    c1: Fp(0x3d17086a9fb936fc),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xdd3b759d0c5b0fa),
                    c1: Fp(0x2b0e767dd35bdfda),
                },
                c1: Fp2 {
                    c0: Fp(0x38af1d5dc3ac6f3e),
                    c1: Fp(0x80d10fafedbe98f),
                },
                c2: Fp2 {
                    c0: Fp(0xad1dd4eb083da47),
                    c1: Fp(0x1ca8311dff4ce349),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x19ff2c168042eb01),
                    c1: Fp(0x4007586da2894e46),
                },
                c1: Fp2 {
                    c0: Fp(0x5aaaf9e3567510a),
                    c1: Fp(0x1a6e418ad84c4619),
                },
                c2: Fp2 {
                    c0: Fp(0x10f896cc9eb1055a),
                    c1: Fp(0x4320112227f80d4),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3bc8fd80f65d02c9),
                    c1: Fp(0x2df2167aed52ddcd),
                },
                c1: Fp2 {
                    c0: Fp(0x1e2041a6d8b5dbbf),
                    c1: Fp(0xc0b9c7468240108),
                },
                c2: Fp2 {
                    c0: Fp(0x7bbabd750c8db9d),
                    c1: Fp(0x1e91996097c9a6d9),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1953c238cdc172e5),
                    c1: Fp(0xe3b4ab21c933459),
                },
                c1: Fp2 {
                    c0: Fp(0x29818ba8470c4fe5),
                    c1: Fp(0x1ae0e695c0f6682e),
                },
                c2: Fp2 {
                    c0: Fp(0x4f2f742ece341b9),
                    c1: Fp(0x361337895c5705bd),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xa1f60ec5d184af1),
                    c1: Fp(0x3a68f25f218c8d02),
                },
                c1: Fp2 {
                    c0: Fp(0x24827bb71a123064),
                    c1: Fp(0x596aa2bee9cc656),
                },
                c2: Fp2 {
                    c0: Fp(0x22247624f7c3ff2c),
                    c1: Fp(0x3794286303ffa006),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3b78586150b23b5d),
                    c1: Fp(0xfab198f88db3e38),
                },
                c1: Fp2 {
                    c0: Fp(0x1405e01d97b2ef1f),
                    c1: Fp(0x3d16c316e3e450e),
                },
                c2: Fp2 {
                    c0: Fp(0x219b40d84bae85f8),
                    c1: Fp(0x18500c542a1eec8a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x493eed46a3a9cfa),
                    c1: Fp(0xd61e31ec6554390),
                },
                c1: Fp2 {
                    c0: Fp(0xfeef9277cf8683a),
                    c1: Fp(0x91fdd8fb889acd0),
                },
                c2: Fp2 {
                    c0: Fp(0x372d0d5088118a97),
                    c1: Fp(0x3d99bede521cd17a),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1987f5100ca21049),
                    c1: Fp(0x1414701a4eb51be8),
                },
                c1: Fp2 {
                    c0: Fp(0x39b1e761da2a3594),
                    c1: Fp(0x408317c7ec7ada0f),
                },
                c2: Fp2 {
                    c0: Fp(0x1ab9dd39f89cc843),
                    c1: Fp(0x3249f35c4bf35541),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3d0558c886ed752c),
                    c1: Fp(0x1d050504926903cc),
                },
                c1: Fp2 {
                    c0: Fp(0x1a74a8847d5a18cc),
                    c1: Fp(0x35b22e1158fbe401),
                },
                c2: Fp2 {
                    c0: Fp(0x197333f74127bb95),
                    c1: Fp(0x6eae05ba4ecf9fa),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xcefb707e333ec0e),
                    c1: Fp(0x3f2cc5547fe4218),
                },
                c1: Fp2 {
                    c0: Fp(0x1c6f4bfe67b4ab84),
                    c1: Fp(0x2e340ebfef02091c),
                },
                c2: Fp2 {
                    c0: Fp(0x13fb7a713c91ba84),
                    c1: Fp(0x3cc0ccb178701160),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x50c7c18fc8a4b98),
                    c1: Fp(0x1a5b07fc191b915c),
                },
                c1: Fp2 {
                    c0: Fp(0x34d906d084fa0611),
                    c1: Fp(0xcbbc8f54de73743),
                },
                c2: Fp2 {
                    c0: Fp(0x25e622bb5d06b99f),
                    c1: Fp(0x23a43c39cbe150),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x42beb89da355cf3),
                    c1: Fp(0x3e71a49e61b9c508),
                },
                c1: Fp2 {
                    c0: Fp(0x21e49e7c4ffaabf1),
                    c1: Fp(0x102ffcbd8f79d5d0),
                },
                c2: Fp2 {
                    c0: Fp(0x24e0406aefcf5f75),
                    c1: Fp(0x3f24e32e613fbc40),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ba5d4b085c0a22c),
                    c1: Fp(0x408c5b24baa5ac60),
                },
                c1: Fp2 {
                    c0: Fp(0x299fae28f95c8c78),
                    c1: Fp(0x39a5e846ce20ebbb),
                },
                c2: Fp2 {
                    c0: Fp(0x20138e2abf79ea0c),
                    c1: Fp(0xe0760764d5aa876),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1c02a1f9f0fc3ac8),
                    c1: Fp(0x3228fae6b9f7b8e4),
                },
                c1: Fp2 {
                    c0: Fp(0x2d4395326f9a1f20),
                    c1: Fp(0x3df313acd576f8d1),
                },
                c2: Fp2 {
                    c0: Fp(0x10be2c1b38a9bbc3),
                    c1: Fp(0x38bc94226c03f0d5),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x366c622262b019bc),
                    c1: Fp(0x7ea4e41d44b7942),
                },
                c1: Fp2 {
                    c0: Fp(0x2ca051041302bc2a),
                    c1: Fp(0x257cf4cecd54b182),
                },
                c2: Fp2 {
                    c0: Fp(0x1510208930eec96a),
                    c1: Fp(0x1ded4a3ea90c518c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x296c27f9aa20efe8),
                    c1: Fp(0x3f273ee4e1c32e64),
                },
                c1: Fp2 {
                    c0: Fp(0x378fc1e844acf273),
                    c1: Fp(0x1d7031a196fcf80a),
                },
                c2: Fp2 {
                    c0: Fp(0xa7cb63c675efcc0),
                    c1: Fp(0x1dba7acbffa128e8),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x39e47df8c163a5fa),
                    c1: Fp(0x3cdcfff6f2d3d363),
                },
                c1: Fp2 {
                    c0: Fp(0x3abad77d5244dc3),
                    c1: Fp(0x1cd3a414c5c9b61d),
                },
                c2: Fp2 {
                    c0: Fp(0x1a5e0474f223aafe),
                    c1: Fp(0x109f8841e31ef65c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x6c522cb19e197cb),
                    c1: Fp(0x5f9e355e8ca9063),
                },
                c1: Fp2 {
                    c0: Fp(0x2f8f2489ff141721),
                    c1: Fp(0x304c95d7f5064cf9),
                },
                c2: Fp2 {
                    c0: Fp(0x3eb9101e7b4829d4),
                    c1: Fp(0x1c8ea4b171932d1d),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1455f2236cde6719),
                    c1: Fp(0x339f4e56ca33e8da),
                },
                c1: Fp2 {
                    c0: Fp(0x40561d93f24edb17),
                    c1: Fp(0x334af21832239bf8),
                },
                c2: Fp2 {
                    c0: Fp(0x415ba34f30213405),
                    c1: Fp(0x34648849e7abce7f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x198f48a59d3f426c),
                    c1: Fp(0x2fa8d6786379a72a),
                },
                c1: Fp2 {
                    c0: Fp(0x2af5b26b7500195),
                    c1: Fp(0x40e59b14fc544dd6),
                },
                c2: Fp2 {
                    c0: Fp(0x387fa190e7006f75),
                    c1: Fp(0x11fdc58256d57edb),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x62949d2a34776c3),
                    c1: Fp(0x4aa24d68132850a),
                },
                c1: Fp2 {
                    c0: Fp(0x377baaa05bd960bf),
                    c1: Fp(0x2df02317deaf3796),
                },
                c2: Fp2 {
                    c0: Fp(0x38de9a3c714487c1),
                    c1: Fp(0x41205958fa142730),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ec319549c034eaa),
                    c1: Fp(0xb77986fbfa5f917),
                },
                c1: Fp2 {
                    c0: Fp(0x39c08f8ef844161a),
                    c1: Fp(0xa9c8d36306b785c),
                },
                c2: Fp2 {
                    c0: Fp(0x20a2a26b50c1e794),
                    c1: Fp(0x60d61d0cbf51d59),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3875bb1287b1cb5e),
                    c1: Fp(0xb26d2ca6185efcb),
                },
                c1: Fp2 {
                    c0: Fp(0xe15851a451c99fb),
                    c1: Fp(0x3db481aeeda8634f),
                },
                c2: Fp2 {
                    c0: Fp(0x4068e59d58229e1d),
                    c1: Fp(0x1eb684543d7f96e7),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1452bf6e8d6b4aed),
                    c1: Fp(0x36fbbfb4184e77c9),
                },
                c1: Fp2 {
                    c0: Fp(0x2a7c077ab94ea3c4),
                    c1: Fp(0x8a410b6ed6e8e54),
                },
                c2: Fp2 {
                    c0: Fp(0xf6518adcfc85764),
                    c1: Fp(0x164703370651fb02),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe8dc17ecc51b3ff),
                    c1: Fp(0x168750ea8f87947d),
                },
                c1: Fp2 {
                    c0: Fp(0xa264abcb9b45815),
                    c1: Fp(0x9284effa33b4cd6),
                },
                c2: Fp2 {
                    c0: Fp(0x63788ff5954da1a),
                    c1: Fp(0x3b52aeaa2b3229d4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x334ba3f21a794059),
                    c1: Fp(0x1c0eece248ce612c),
                },
                c1: Fp2 {
                    c0: Fp(0x33d44d7f27b7353),
                    c1: Fp(0x134742dd3a7f4708),
                },
                c2: Fp2 {
                    c0: Fp(0x3519887cf20ba136),
                    c1: Fp(0x2d504b0be7467fac),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c4e1a9fd9ce5b45),
                    c1: Fp(0x1cfe2273186a6fa2),
                },
                c1: Fp2 {
                    c0: Fp(0x2d14a5af5bfc703b),
                    c1: Fp(0x30bbb02a767aa27a),
                },
                c2: Fp2 {
                    c0: Fp(0x413093e31166589a),
                    c1: Fp(0x13638e116fa7a7ed),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x11defcf5e620a367),
                    c1: Fp(0x329b8ad75dbc1353),
                },
                c1: Fp2 {
                    c0: Fp(0x1a05086a84346df0),
                    c1: Fp(0x17633194ce79339b),
                },
                c2: Fp2 {
                    c0: Fp(0x299b7d87932c2b2),
                    c1: Fp(0x408fd90d82b3ae0b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x21123a7d59487eac),
                    c1: Fp(0x160576a9a57356b4),
                },
                c1: Fp2 {
                    c0: Fp(0x376b0cbab8c9c0fc),
                    c1: Fp(0x2cd6dac4c4f1e78a),
                },
                c2: Fp2 {
                    c0: Fp(0xcad8e636746ffb),
                    c1: Fp(0x204f95d6d19931be),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x39d4b80e09558d98),
                    c1: Fp(0x2cc395633db09334),
                },
                c1: Fp2 {
                    c0: Fp(0x30f4b4d74e548c6c),
                    c1: Fp(0x3e0fe545ff5f05b2),
                },
                c2: Fp2 {
                    c0: Fp(0xd092ba687546031),
                    c1: Fp(0x15acc01c7f95e84),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x383f2d2fc5a81ef7),
                    c1: Fp(0x13d46c4420605c3f),
                },
                c1: Fp2 {
                    c0: Fp(0x6cd3600f1043394),
                    c1: Fp(0x3451efd53506794f),
                },
                c2: Fp2 {
                    c0: Fp(0x1291aabbd3108c9d),
                    c1: Fp(0x773818506e778fe),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1f3888c2a72c36bd),
                    c1: Fp(0x22db7b2cdbea5dd1),
                },
                c1: Fp2 {
                    c0: Fp(0x1ca3ff041551298a),
                    c1: Fp(0x3108c84576c59627),
                },
                c2: Fp2 {
                    c0: Fp(0x34bf2d868c48c27e),
                    c1: Fp(0x121a5b715b55003c),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28bf3c5af38d9554),
                    c1: Fp(0x357c62d7dbec3860),
                },
                c1: Fp2 {
                    c0: Fp(0x1d2ddb20baa72438),
                    c1: Fp(0x52620489828fb50),
                },
                c2: Fp2 {
                    c0: Fp(0x2af2649ecd6f6056),
                    c1: Fp(0x37ba29549030353b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x364ef466efde59d9),
                    c1: Fp(0x11db385be6438525),
                },
                c1: Fp2 {
                    c0: Fp(0x1f53e297811616fc),
                    c1: Fp(0x4134eee5d588675),
                },
                c2: Fp2 {
                    c0: Fp(0x568ce4c8fdfa320),
                    c1: Fp(0x3e5ac8853c9003a7),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2a47f6a8203f312b),
                    c1: Fp(0x28f97b42161241f7),
                },
                c1: Fp2 {
                    c0: Fp(0x9c9d12bf3705f39),
                    c1: Fp(0x133896a48028c8),
                },
                c2: Fp2 {
                    c0: Fp(0x207a6f11ca2cd07),
                    c1: Fp(0xa73633e24195fcc),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x76f0ca8ea4b2bc2),
                    c1: Fp(0x32b1f972cf84cdd2),
                },
                c1: Fp2 {
                    c0: Fp(0x3c1c173787ffa305),
                    c1: Fp(0x331269da4397127),
                },
                c2: Fp2 {
                    c0: Fp(0x3347bd12bda1c9e0),
                    c1: Fp(0x26b59b2cf11a78ea),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2d407481b6dbfc34),
                    c1: Fp(0x37db092246d84ed1),
                },
                c1: Fp2 {
                    c0: Fp(0x1165db712def3cc),
                    c1: Fp(0x401315edf5a7f46),
                },
                c2: Fp2 {
                    c0: Fp(0x167adc94f5e0e758),
                    c1: Fp(0x486ca8d515c22b6),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x262053d682dc5c58),
                    c1: Fp(0x1241032dcf752cba),
                },
                c1: Fp2 {
                    c0: Fp(0x370e5ef84fbba40e),
                    c1: Fp(0x1e7c1ce3ecfb732a),
                },
                c2: Fp2 {
                    c0: Fp(0x4c861b99d6fc001),
                    c1: Fp(0x1caff9feea5b17da),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x387aa4b688eef66),
                    c1: Fp(0x200e321306df883e),
                },
                c1: Fp2 {
                    c0: Fp(0x28ed795e110f0da6),
                    c1: Fp(0x5ab6a43acd1a86),
                },
                c2: Fp2 {
                    c0: Fp(0x3eb83dd0192ba75d),
                    c1: Fp(0x3f715c88ca0969e1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xcacd3e51dbb8f55),
                    c1: Fp(0x30d3d0a955c831c8),
                },
                c1: Fp2 {
                    c0: Fp(0x3842953fecf44c6c),
                    c1: Fp(0x8a954dd8111b552),
                },
                c2: Fp2 {
                    c0: Fp(0x25aac2e0c43f5412),
                    c1: Fp(0x2b54a825b0c7504e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x396416f4645ad28d),
                    c1: Fp(0x1dd93281be9aaab2),
                },
                c1: Fp2 {
                    c0: Fp(0x135ca2ace64a43b2),
                    c1: Fp(0x8b7e38f97b9c930),
                },
                c2: Fp2 {
                    c0: Fp(0xd9a02459066b30a),
                    c1: Fp(0x22b0b9b5f33f0256),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x9209a6c4ab962b0),
                    c1: Fp(0x4155c0432771499),
                },
                c1: Fp2 {
                    c0: Fp(0xd199932635e5a93),
                    c1: Fp(0x327a87633f093eda),
                },
                c2: Fp2 {
                    c0: Fp(0x30a7a27900524520),
                    c1: Fp(0xb5b730d34931e2e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x26c4b000498763ba),
                    c1: Fp(0x40680b5926f1bf1d),
                },
                c1: Fp2 {
                    c0: Fp(0x4870c32b57673aa),
                    c1: Fp(0xb605c087809f51a),
                },
                c2: Fp2 {
                    c0: Fp(0x1b9c57ada851202d),
                    c1: Fp(0xdbe58cc32e43f7e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1c00a58a7f51dfa),
                    c1: Fp(0x17ff9b8e616c767d),
                },
                c1: Fp2 {
                    c0: Fp(0x26c62aa0e9c2f2da),
                    c1: Fp(0x226914dfa3148b82),
                },
                c2: Fp2 {
                    c0: Fp(0x3a509d1cc1d228ca),
                    c1: Fp(0xcf01b9d4ada80cb),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a16421487b2d139),
                    c1: Fp(0x2ef8ab0b3eabe36d),
                },
                c1: Fp2 {
                    c0: Fp(0x106cd3bd891e5ca),
                    c1: Fp(0x2da3b5c2f7b5f731),
                },
                c2: Fp2 {
                    c0: Fp(0x326e1975c861dc82),
                    c1: Fp(0x14577b9fe9865dae),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25a1e02159c57b95),
                    c1: Fp(0x3fa4f09b5f001bd7),
                },
                c1: Fp2 {
                    c0: Fp(0x6d3231add627512),
                    c1: Fp(0xf8937ad99c78a27),
                },
                c2: Fp2 {
                    c0: Fp(0x2c6499fad1c83f8d),
                    c1: Fp(0x358fda2c2212c3fb),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28bfc1ac20542a0d),
                    c1: Fp(0x1ca9956596f1271c),
                },
                c1: Fp2 {
                    c0: Fp(0x2fc47aaffd0e8523),
                    c1: Fp(0x5fd5df6268a4a2d),
                },
                c2: Fp2 {
                    c0: Fp(0x29f57b5d038769f8),
                    c1: Fp(0x3e03b243837ce742),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x5aabe88119b2fff),
                    c1: Fp(0x1d921384a65a80e7),
                },
                c1: Fp2 {
                    c0: Fp(0x5f578be337f9b2a),
                    c1: Fp(0x2ddd80e3aa58a616),
                },
                c2: Fp2 {
                    c0: Fp(0xfc40e02c1fd4678),
                    c1: Fp(0x247a3cfa61950ca7),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2d7d3a236eec9e09),
                    c1: Fp(0x3bd3e9b4f0fdb766),
                },
                c1: Fp2 {
                    c0: Fp(0x3cab56c9582949a1),
                    c1: Fp(0x3d9999b576dd2868),
                },
                c2: Fp2 {
                    c0: Fp(0x3ddfcf5895c0ba6b),
                    c1: Fp(0xaa063d770599979),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xee23fb367952402),
                    c1: Fp(0x15aa622a9653c190),
                },
                c1: Fp2 {
                    c0: Fp(0x1d5576ed8064a097),
                    c1: Fp(0x2274a35b594324ce),
                },
                c2: Fp2 {
                    c0: Fp(0x26cab7503b2f717),
                    c1: Fp(0x381d7419e5045ec2),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xb5d23c085717d8c),
                    c1: Fp(0xbf0ff95d6742e4c),
                },
                c1: Fp2 {
                    c0: Fp(0x3aa2d0f787e660d2),
                    c1: Fp(0x965e34d612ab641),
                },
                c2: Fp2 {
                    c0: Fp(0x3609f4006a80f76),
                    c1: Fp(0xec5e1ea5da2a92),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x33005078ffc06086),
                    c1: Fp(0x1dcd337f3eeeebbc),
                },
                c1: Fp2 {
                    c0: Fp(0x1f2d0d2a7b9db5f5),
                    c1: Fp(0x40309f435ba7ddd9),
                },
                c2: Fp2 {
                    c0: Fp(0x669e1d4d783e8de),
                    c1: Fp(0x1e1ba2081d4f520f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x32a22a2c9ab7ed1a),
                    c1: Fp(0x2dc1cbbe2fe194f8),
                },
                c1: Fp2 {
                    c0: Fp(0x36a74abd154d8424),
                    c1: Fp(0x7e08ec84448d367),
                },
                c2: Fp2 {
                    c0: Fp(0x3ab07c33bb61bb5f),
                    c1: Fp(0x1186ed7e9cd7e85),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2f3449b7f597f00a),
                    c1: Fp(0x1f66b45b2b6f87d0),
                },
                c1: Fp2 {
                    c0: Fp(0x294af011906ca60d),
                    c1: Fp(0x337fa73c07ba5181),
                },
                c2: Fp2 {
                    c0: Fp(0xf3b0caffba52360),
                    c1: Fp(0x1beded46b28bdcef),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3bb22a481d0cb232),
                    c1: Fp(0x7997208ff246dad),
                },
                c1: Fp2 {
                    c0: Fp(0xac32e6a226baf97),
                    c1: Fp(0x11cebad135fb0edf),
                },
                c2: Fp2 {
                    c0: Fp(0x192fa3ade10484ab),
                    c1: Fp(0x3bf9fe8080f27052),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x193c522346ea29a1),
                    c1: Fp(0xae1ba1bf7a30c73),
                },
                c1: Fp2 {
                    c0: Fp(0x1d358b8dba3734ee),
                    c1: Fp(0x1465544507797a6b),
                },
                c2: Fp2 {
                    c0: Fp(0x1b614e1f62dada0f),
                    c1: Fp(0x23974cf061243099),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x20729a26bc46b59b),
                    c1: Fp(0x455aa56aa31876d),
                },
                c1: Fp2 {
                    c0: Fp(0x8d857c28206e2fd),
                    c1: Fp(0x3b4c041c40fb5c53),
                },
                c2: Fp2 {
                    c0: Fp(0x197a33627d4bc790),
                    c1: Fp(0x2b91173a628f70e8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1122bc600fd6d90b),
                    c1: Fp(0x274e4e0f6031d142),
                },
                c1: Fp2 {
                    c0: Fp(0x25b11fca35ebb88b),
                    c1: Fp(0x7287f6ccae87fa9),
                },
                c2: Fp2 {
                    c0: Fp(0x203b736416b61874),
                    c1: Fp(0x534a9f194e7743d),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ca8119c8b9e59b4),
                    c1: Fp(0xe57843d39a39916),
                },
                c1: Fp2 {
                    c0: Fp(0xfd64b2937ae12ac),
                    c1: Fp(0x739435125544deb),
                },
                c2: Fp2 {
                    c0: Fp(0x1af5362430a05b95),
                    c1: Fp(0x3a50ad1c1ec6a7af),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x40a2eb146a879993),
                    c1: Fp(0x33666791836dabfd),
                },
                c1: Fp2 {
                    c0: Fp(0x17cd1e5b67e11eea),
                    c1: Fp(0x1e3703f13f2911e),
                },
                c2: Fp2 {
                    c0: Fp(0x200f9993316edae8),
                    c1: Fp(0x129f3fffada044a0),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28ca5134e159e9f5),
                    c1: Fp(0x74be3311bd00830),
                },
                c1: Fp2 {
                    c0: Fp(0x845fd3533319c50),
                    c1: Fp(0x3ad0f64d2e3ca313),
                },
                c2: Fp2 {
                    c0: Fp(0x123407c23a623e48),
                    c1: Fp(0x210bd09864cc9643),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3b1f4fc97805c47c),
                    c1: Fp(0xe89c7db17217bf9),
                },
                c1: Fp2 {
                    c0: Fp(0x9f4fb50a008d3da),
                    c1: Fp(0x2cb06800273255dd),
                },
                c2: Fp2 {
                    c0: Fp(0x17db8f55c766c2f0),
                    c1: Fp(0x2c0972cafb631120),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x120ae61b5f064599),
                    c1: Fp(0x35747ad440ffd73c),
                },
                c1: Fp2 {
                    c0: Fp(0x125741c1c2c8d519),
                    c1: Fp(0x1b4b88d377c8297a),
                },
                c2: Fp2 {
                    c0: Fp(0xfd2af5afb2d75a5),
                    c1: Fp(0x1cd0762d664a18b5),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x4bd31d05c26f084),
                    c1: Fp(0x3b33aeda7eddfa6),
                },
                c1: Fp2 {
                    c0: Fp(0x2c79e9ca8fe3e724),
                    c1: Fp(0x33141618820272e7),
                },
                c2: Fp2 {
                    c0: Fp(0x372a9be2e477de9f),
                    c1: Fp(0x1812ae657a8c3ba6),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2e25a5686455ebe9),
                    c1: Fp(0x178d85ae8aad6c48),
                },
                c1: Fp2 {
                    c0: Fp(0x3ed5cdb4293f4fe8),
                    c1: Fp(0xcc34cc68fc26b30),
                },
                c2: Fp2 {
                    c0: Fp(0x2cf56e2ebc6ac70c),
                    c1: Fp(0x82b87316722253c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x340a8f2915b405),
                    c1: Fp(0x327a53a9035f23d3),
                },
                c1: Fp2 {
                    c0: Fp(0x2c4db617cb548426),
                    c1: Fp(0x1f6a3b3a5ac24986),
                },
                c2: Fp2 {
                    c0: Fp(0x224708b0173b6558),
                    c1: Fp(0x2741f9dcb2a0493b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x178f9b72b35f741b),
                    c1: Fp(0x21ee87cbb5d42f35),
                },
                c1: Fp2 {
                    c0: Fp(0x3b170908e30e88ba),
                    c1: Fp(0x2f2d34a5596e8c64),
                },
                c2: Fp2 {
                    c0: Fp(0x398eacc60eba4217),
                    c1: Fp(0x26d330879fed54bf),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28ee45f5cc2fa3b),
                    c1: Fp(0x19a30980cbae7329),
                },
                c1: Fp2 {
                    c0: Fp(0x1d3ab896e735e431),
                    c1: Fp(0x197f2ea9ef1b7db6),
                },
                c2: Fp2 {
                    c0: Fp(0x185ec3a2a6244b5b),
                    c1: Fp(0x2b744ec39ccf9336),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x35927217c7e6bdc0),
                    c1: Fp(0x286919f8205cf176),
                },
                c1: Fp2 {
                    c0: Fp(0x161b36110ba184e0),
                    c1: Fp(0x2c5be9901e61ba36),
                },
                c2: Fp2 {
                    c0: Fp(0x23d5eaae31a69bca),
                    c1: Fp(0x7c798a6661ed1a3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1c08e17cd5a374fa),
                    c1: Fp(0x3bd8c74aa0e30403),
                },
                c1: Fp2 {
                    c0: Fp(0x1602a3378cded6fc),
                    c1: Fp(0x415c4977a7e2ea82),
                },
                c2: Fp2 {
                    c0: Fp(0x1146642656a0c4e7),
                    c1: Fp(0x25ebaddd18aff04e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1aa63f91b9874aee),
                    c1: Fp(0x2cd4fefe040d307),
                },
                c1: Fp2 {
                    c0: Fp(0x3894c30bac83de39),
                    c1: Fp(0x1c07d99063c1d254),
                },
                c2: Fp2 {
                    c0: Fp(0x2feb998f53056a6e),
                    c1: Fp(0x3efc9ad2945450bb),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3417f128c1c9c766),
                    c1: Fp(0x1aa3feb86c162726),
                },
                c1: Fp2 {
                    c0: Fp(0x9e4386a6610eef2),
                    c1: Fp(0xb034ecc21a06c90),
                },
                c2: Fp2 {
                    c0: Fp(0x24a6425ac5ef2727),
                    c1: Fp(0x2d4453f3124859c5),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x330a3e5c362d22ea),
                    c1: Fp(0x160ba1da33cdb3b9),
                },
                c1: Fp2 {
                    c0: Fp(0x2e6e003b05017606),
                    c1: Fp(0x298f42e7d0c975d9),
                },
                c2: Fp2 {
                    c0: Fp(0x1b43c7edf5d929e),
                    c1: Fp(0x2b010fb77f4188bc),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x315d0615df2843c3),
                    c1: Fp(0x2d119a71048396b1),
                },
                c1: Fp2 {
                    c0: Fp(0x3d7d038f04dd6147),
                    c1: Fp(0x1648abb691d404ca),
                },
                c2: Fp2 {
                    c0: Fp(0xb527edae3a413c5),
                    c1: Fp(0x1568c9cf7f39dbad),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2ae33929d9f5f1fb),
                    c1: Fp(0x47de99ee98fd139),
                },
                c1: Fp2 {
                    c0: Fp(0x2fb14a85ca798111),
                    c1: Fp(0x3cc8fa121d313889),
                },
                c2: Fp2 {
                    c0: Fp(0x24ca5d1cef51a28c),
                    c1: Fp(0x33105ec92e95ffb1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2cb3ddb2d67dd19a),
                    c1: Fp(0x89491b07772b0c6),
                },
                c1: Fp2 {
                    c0: Fp(0x33e0d76b2598f0ec),
                    c1: Fp(0x1935d79068ebf641),
                },
                c2: Fp2 {
                    c0: Fp(0x1369e91c010a4606),
                    c1: Fp(0x27b77fc4f33956ff),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3b91e98810aa02a4),
                    c1: Fp(0x2284ab07def5958),
                },
                c1: Fp2 {
                    c0: Fp(0xa5fb16df921230f),
                    c1: Fp(0xb17f493346c76ca),
                },
                c2: Fp2 {
                    c0: Fp(0xbf06dadf7c4667e),
                    c1: Fp(0x3445b83e6436037a),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x23da3ae7468eb330),
                    c1: Fp(0x21b90ea00c94bbe9),
                },
                c1: Fp2 {
                    c0: Fp(0x295b9cfac385f0e0),
                    c1: Fp(0x7f14a71d91038d4),
                },
                c2: Fp2 {
                    c0: Fp(0x20e3702b4aed5978),
                    c1: Fp(0x3619cf4f04431a0f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xff427424e34b76c),
                    c1: Fp(0x2dd6d1bfba0b936e),
                },
                c1: Fp2 {
                    c0: Fp(0x2eeea57afa87c368),
                    c1: Fp(0x65131d89539fc8b),
                },
                c2: Fp2 {
                    c0: Fp(0x3ff915fc078f2ab7),
                    c1: Fp(0xb9ad96795b097cb),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x6acfbdb2bf99939),
                    c1: Fp(0x184e46446c80b59c),
                },
                c1: Fp2 {
                    c0: Fp(0x342fb947fd106258),
                    c1: Fp(0xf99f46374dc1aab),
                },
                c2: Fp2 {
                    c0: Fp(0x14efc5499dbe58e3),
                    c1: Fp(0x15b3719661c71bfc),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x616fe4a1676b585),
                    c1: Fp(0x2bcb3a47453a6b5f),
                },
                c1: Fp2 {
                    c0: Fp(0x26de220e617bf4f5),
                    c1: Fp(0x16817dadb5cb122d),
                },
                c2: Fp2 {
                    c0: Fp(0xc959a2ce05a7380),
                    c1: Fp(0x62e40073c819722),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x11dcf0309c1b7cd),
                    c1: Fp(0x328e291e678fb4eb),
                },
                c1: Fp2 {
                    c0: Fp(0x30cf946b2ee96c6f),
                    c1: Fp(0x2a57cfc1e95a0b0f),
                },
                c2: Fp2 {
                    c0: Fp(0x137556b9c39a2bae),
                    c1: Fp(0x349d57e0aad6af76),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xc33d85ce9bf7cb9),
                    c1: Fp(0x110bda24af2da9e5),
                },
                c1: Fp2 {
                    c0: Fp(0x30827bc30acb1065),
                    c1: Fp(0x2c5a3a5c7c00ec30),
                },
                c2: Fp2 {
                    c0: Fp(0x3dabbd4e5eda60ae),
                    c1: Fp(0x2a370de12d9c0255),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2210c81ea8024428),
                    c1: Fp(0x3958eda28bdbf747),
                },
                c1: Fp2 {
                    c0: Fp(0x2a213ce5354799f6),
                    c1: Fp(0x1910e7bd5940df56),
                },
                c2: Fp2 {
                    c0: Fp(0xb6bc64f068a99b0),
                    c1: Fp(0x18afbb9fa62bd98c),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe5ba39db0e98b4b),
                    c1: Fp(0x3861ccd659c82757),
                },
                c1: Fp2 {
                    c0: Fp(0x2d3fc54b93742e79),
                    c1: Fp(0x5d5dca65adf2e59),
                },
                c2: Fp2 {
                    c0: Fp(0x1e031e806db22671),
                    c1: Fp(0x197ad5849e6b4437),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x10ab026ed43a806),
                    c1: Fp(0x403b7798fd0b5dad),
                },
                c1: Fp2 {
                    c0: Fp(0x1ec2a836b6839b8b),
                    c1: Fp(0x3ee7a37b518e44b9),
                },
                c2: Fp2 {
                    c0: Fp(0x3a57764aa7655402),
                    c1: Fp(0x414c1d66a3fde26a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1dded56454d72201),
                    c1: Fp(0x2efe5d8bff2f17e5),
                },
                c1: Fp2 {
                    c0: Fp(0xcbb4904558d671c),
                    c1: Fp(0x32a959f080232774),
                },
                c2: Fp2 {
                    c0: Fp(0x1a0a4d995af3963b),
                    c1: Fp(0x233c79c25aa542ac),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x16e08db4784643e7),
                    c1: Fp(0x2c6e1f5a22991624),
                },
                c1: Fp2 {
                    c0: Fp(0x31441f65dc78a2e0),
                    c1: Fp(0x14ef4b0e2275ec7e),
                },
                c2: Fp2 {
                    c0: Fp(0x14ae80e16badfcd),
                    c1: Fp(0x71e26e3c1654177),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1b9fbd204661a740),
                    c1: Fp(0x3c6ac087336b463d),
                },
                c1: Fp2 {
                    c0: Fp(0x21338bc589430860),
                    c1: Fp(0x125b8696716031be),
                },
                c2: Fp2 {
                    c0: Fp(0x1e31604c17543a13),
                    c1: Fp(0x2d1f3eaf59029174),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x166dbe32062af235),
                    c1: Fp(0x3e9712a462b7e355),
                },
                c1: Fp2 {
                    c0: Fp(0x2ff18d83d56fccdb),
                    c1: Fp(0x36347f493c68d555),
                },
                c2: Fp2 {
                    c0: Fp(0x20cfd9c418e09ff4),
                    c1: Fp(0x3b096bf392167505),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x289b972fe72a6ed7),
                    c1: Fp(0x345588e453ee5b37),
                },
                c1: Fp2 {
                    c0: Fp(0x2b8354c6b463009d),
                    c1: Fp(0x357585fe9ec95032),
                },
                c2: Fp2 {
                    c0: Fp(0x2a7f2577bbbb2805),
                    c1: Fp(0x82af3e284f4096f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2ef3e59ac966fdda),
                    c1: Fp(0x23346b87f2aa42f9),
                },
                c1: Fp2 {
                    c0: Fp(0x36f02985df99f088),
                    c1: Fp(0x1951b16bdf91c406),
                },
                c2: Fp2 {
                    c0: Fp(0x233139c4fb3138ed),
                    c1: Fp(0x338c5e0bbdb4a6dd),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a80b9f4f1688b2b),
                    c1: Fp(0x33b1f2c641cc0f36),
                },
                c1: Fp2 {
                    c0: Fp(0x215aafaee2b499b8),
                    c1: Fp(0x3371b6a51025e336),
                },
                c2: Fp2 {
                    c0: Fp(0x17575823d4d63fc7),
                    c1: Fp(0xf75541b866811b6),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x517eb75c33ce322),
                    c1: Fp(0x2b841f189d48077),
                },
                c1: Fp2 {
                    c0: Fp(0x9a32111f720faa9),
                    c1: Fp(0x181a9fd6f30f4ab0),
                },
                c2: Fp2 {
                    c0: Fp(0x2a3930e02285ed2d),
                    c1: Fp(0x127c45c8aa3c8b7f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x24e4e586cede78e1),
                    c1: Fp(0xc8bd1e3a2823b2b),
                },
                c1: Fp2 {
                    c0: Fp(0x3cf252c9addaf9eb),
                    c1: Fp(0x3a93417bee4bac31),
                },
                c2: Fp2 {
                    c0: Fp(0x39731e52d147d388),
                    c1: Fp(0x2ea2f0dc67f7278c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28ec8e196989494a),
                    c1: Fp(0x417c617a06cc07d6),
                },
                c1: Fp2 {
                    c0: Fp(0x37fc16148b99a681),
                    c1: Fp(0x29b403c0b20e6382),
                },
                c2: Fp2 {
                    c0: Fp(0x1b2aae027973ef8),
                    c1: Fp(0x3f0da0cdf63e7da9),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3dc878986f31d250),
                    c1: Fp(0x28f082cba82e2d72),
                },
                c1: Fp2 {
                    c0: Fp(0x37888c2c02cdf81e),
                    c1: Fp(0x2e255f20169390d),
                },
                c2: Fp2 {
                    c0: Fp(0x3a48e5998cb0d0be),
                    c1: Fp(0x1272e76540784592),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xfac95f4a55302f1),
                    c1: Fp(0x1e23f699dfe2a8a8),
                },
                c1: Fp2 {
                    c0: Fp(0x347569c9a220786f),
                    c1: Fp(0x39bb6b4a07288df6),
                },
                c2: Fp2 {
                    c0: Fp(0x3c27b243c4b400f2),
                    c1: Fp(0x66c5b937123fe78),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c85433c120a5667),
                    c1: Fp(0x1827cfebeab3b4ba),
                },
                c1: Fp2 {
                    c0: Fp(0x29a47d124db9588d),
                    c1: Fp(0x2f9143f11435982),
                },
                c2: Fp2 {
                    c0: Fp(0x2a1c55b8269f631c),
                    c1: Fp(0x347c6c40a81c9d60),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ad260b3d48efaff),
                    c1: Fp(0xe5901d22d592958),
                },
                c1: Fp2 {
                    c0: Fp(0x9fb17591497537f),
                    c1: Fp(0x33efa29a96967ac2),
                },
                c2: Fp2 {
                    c0: Fp(0x14ce7a7564531f46),
                    c1: Fp(0x3ee64ca3c7be2ffb),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3100cba2dc9bf24b),
                    c1: Fp(0x38cd066fdff046b0),
                },
                c1: Fp2 {
                    c0: Fp(0x239828c8ade271dc),
                    c1: Fp(0x1c70a2bce143a340),
                },
                c2: Fp2 {
                    c0: Fp(0x19113f6930720604),
                    c1: Fp(0xb5e3dff5492f3ee),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1eb4098200c057e7),
                    c1: Fp(0x371df4dea153fa6d),
                },
                c1: Fp2 {
                    c0: Fp(0x179e84a749c91888),
                    c1: Fp(0xc202097ee1c75e),
                },
                c2: Fp2 {
                    c0: Fp(0x27fd28f735a7c04f),
                    c1: Fp(0x26ec1a6c1cef46e8),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x161ca90789efa7c4),
                    c1: Fp(0x3e80b43a3c618b10),
                },
                c1: Fp2 {
                    c0: Fp(0x40babc783012fbc5),
                    c1: Fp(0x2aaf75a5edd3f50c),
                },
                c2: Fp2 {
                    c0: Fp(0x119b441b5a980d1e),
                    c1: Fp(0x3427f3472223763e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2594d84e419d9cd1),
                    c1: Fp(0xf0de490f459e04a),
                },
                c1: Fp2 {
                    c0: Fp(0x1343ac16ed8f9244),
                    c1: Fp(0x19e9e8318a2f98f1),
                },
                c2: Fp2 {
                    c0: Fp(0x2f638f1ff29c79aa),
                    c1: Fp(0x3761f2f007cb4c90),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3edf535a4a880755),
                    c1: Fp(0x151c4799be06665b),
                },
                c1: Fp2 {
                    c0: Fp(0x2e23598c85b657a2),
                    c1: Fp(0x9713821c7005b19),
                },
                c2: Fp2 {
                    c0: Fp(0x2ab2533ee5958ef3),
                    c1: Fp(0x2803a718dd566f56),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x13ac54c7823e8830),
                    c1: Fp(0x18c8014a0673c2ae),
                },
                c1: Fp2 {
                    c0: Fp(0x2186ab519ac9af56),
                    c1: Fp(0x16d896a6bbf0b6c0),
                },
                c2: Fp2 {
                    c0: Fp(0x1ea0ac117e3d31dc),
                    c1: Fp(0x20f75f4459f687dc),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x11894a01ba9fc21c),
                    c1: Fp(0x1254f02c99086ccf),
                },
                c1: Fp2 {
                    c0: Fp(0xc3c9ab46dbc5b74),
                    c1: Fp(0x43f6116251be1d),
                },
                c2: Fp2 {
                    c0: Fp(0x19632c705c496d10),
                    c1: Fp(0x22df7ab58493d719),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x27ec803cddb3624c),
                    c1: Fp(0xff16695633007ec),
                },
                c1: Fp2 {
                    c0: Fp(0xfb089c3af60e0b5),
                    c1: Fp(0x40e5d2ee0089fb47),
                },
                c2: Fp2 {
                    c0: Fp(0x25960f04fd153c22),
                    c1: Fp(0xe4e82ed88c82fd6),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2308e28be0dd5655),
                    c1: Fp(0x5166ba0c1261be9),
                },
                c1: Fp2 {
                    c0: Fp(0x3cd97fa93947ac9f),
                    c1: Fp(0x76efa4d4ad85cb3),
                },
                c2: Fp2 {
                    c0: Fp(0xb0cf7a51b46ea06),
                    c1: Fp(0x3f15eaf51c859d07),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x19c899ecfeab0bdf),
                    c1: Fp(0x12c7be2b1fd3e178),
                },
                c1: Fp2 {
                    c0: Fp(0x38d70b0b4e7465be),
                    c1: Fp(0x518b329e3a42b93),
                },
                c2: Fp2 {
                    c0: Fp(0x33ccebd05e0c11dd),
                    c1: Fp(0x3051d4611c01f84e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x37de56da98f63f77),
                    c1: Fp(0x4126ead1eeb04c3e),
                },
                c1: Fp2 {
                    c0: Fp(0x1ce447fb7db5954e),
                    c1: Fp(0x1d7e40e2ffcde330),
                },
                c2: Fp2 {
                    c0: Fp(0x3c459dd5509d0079),
                    c1: Fp(0x2a1760371f06252d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x348f3af9ba96b463),
                    c1: Fp(0xa8af9b161739de8),
                },
                c1: Fp2 {
                    c0: Fp(0x206a1d8427c0614a),
                    c1: Fp(0x2ed656be560fd9f),
                },
                c2: Fp2 {
                    c0: Fp(0x383910f5c7c5aca3),
                    c1: Fp(0x3a0976ba9dc39af7),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x287ac0ade52d111),
                    c1: Fp(0x9b021ec92770310),
                },
                c1: Fp2 {
                    c0: Fp(0x39422eb5d5f01606),
                    c1: Fp(0x3e68be23f5e24fff),
                },
                c2: Fp2 {
                    c0: Fp(0x33cbd6b3a48350d9),
                    c1: Fp(0x29bb1c2cca93252d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x156d04286c6da481),
                    c1: Fp(0x2d7f6cce8a686a56),
                },
                c1: Fp2 {
                    c0: Fp(0x2ffd12de16bd3268),
                    c1: Fp(0xed2e01f3f91d635),
                },
                c2: Fp2 {
                    c0: Fp(0x27c650c0ee11db52),
                    c1: Fp(0x28b828ee2eecce3c),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2d06679cf7aae757),
                    c1: Fp(0xb19b1b24496a233),
                },
                c1: Fp2 {
                    c0: Fp(0x2397efd958d9884d),
                    c1: Fp(0x12de2c08e36e7e7a),
                },
                c2: Fp2 {
                    c0: Fp(0x39ff41337e5fa3d4),
                    c1: Fp(0x11f7efe9edca1a9f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3b409f0e40562a94),
                    c1: Fp(0x6574d1df5ce3bd6),
                },
                c1: Fp2 {
                    c0: Fp(0x20b20c476704195d),
                    c1: Fp(0x2cb79e7a0cd70176),
                },
                c2: Fp2 {
                    c0: Fp(0x25a3092d1b17f2ea),
                    c1: Fp(0xe0a8031abc3f417),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3668f671e4aa2fc8),
                    c1: Fp(0x3f83531bfca542c1),
                },
                c1: Fp2 {
                    c0: Fp(0x22146f3f0c7c17ac),
                    c1: Fp(0x102d1170a701a8d5),
                },
                c2: Fp2 {
                    c0: Fp(0x1493b8c8e5c2e328),
                    c1: Fp(0x3c3e011a76fa4291),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3625af298cff46cf),
                    c1: Fp(0x2f65b09f967599a6),
                },
                c1: Fp2 {
                    c0: Fp(0x3e907d441f3b2515),
                    c1: Fp(0x23460e1fa8098007),
                },
                c2: Fp2 {
                    c0: Fp(0x1b3c76a36803b06e),
                    c1: Fp(0x1be262f1f5e63c26),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25cb4caf1e4fb440),
                    c1: Fp(0x2a6e54d596615b0b),
                },
                c1: Fp2 {
                    c0: Fp(0x3f02a5667054c3e9),
                    c1: Fp(0x1768fa8f9a59df44),
                },
                c2: Fp2 {
                    c0: Fp(0x5f5ec8b893551bd),
                    c1: Fp(0x2780587bac38535),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1aef2e611e7104b0),
                    c1: Fp(0x233cbc500158092a),
                },
                c1: Fp2 {
                    c0: Fp(0x3d39d044dd2e5443),
                    c1: Fp(0x3da3443f054a06f6),
                },
                c2: Fp2 {
                    c0: Fp(0x9ef20d5332e8fde),
                    c1: Fp(0x1cb86b46f12de4f6),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xb6ab9c9bdfd7810),
                    c1: Fp(0x3ef047360af542c2),
                },
                c1: Fp2 {
                    c0: Fp(0x3914215463ec0ca6),
                    c1: Fp(0x2d58c55147a3b78b),
                },
                c2: Fp2 {
                    c0: Fp(0x22fe5ae542edffe0),
                    c1: Fp(0x342b80fc14bd377d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x929afb5b62174cb),
                    c1: Fp(0x2cb96fe613588cb4),
                },
                c1: Fp2 {
                    c0: Fp(0x3fb781e5328b400b),
                    c1: Fp(0x2cde0432027e6745),
                },
                c2: Fp2 {
                    c0: Fp(0x48646c532719d4),
                    c1: Fp(0xad34becd64303dc),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2b1c74b8b01bef12),
                    c1: Fp(0x28d8447f62dc619e),
                },
                c1: Fp2 {
                    c0: Fp(0x3062906311d0e49),
                    c1: Fp(0x3939c583207b3916),
                },
                c2: Fp2 {
                    c0: Fp(0xcdf4f9cf62b8b1),
                    c1: Fp(0xf77b2b38a3aa8e6),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x356eb67cd8a773ea),
                    c1: Fp(0x1c430be72cf24e49),
                },
                c1: Fp2 {
                    c0: Fp(0x6e974cfb65b9d3f),
                    c1: Fp(0x3d55c9789f87841e),
                },
                c2: Fp2 {
                    c0: Fp(0x34a8b1afcb9299e5),
                    c1: Fp(0x119730fc46772b29),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ec3bce93c59e292),
                    c1: Fp(0x1b68c530fb300f45),
                },
                c1: Fp2 {
                    c0: Fp(0x15115eef6e9d093a),
                    c1: Fp(0x1dfbee23696ef2b1),
                },
                c2: Fp2 {
                    c0: Fp(0x11579cfce4998b27),
                    c1: Fp(0x94ec8da20931e6d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x21cfbf7eeb2af312),
                    c1: Fp(0x2bacfa05d69d1f06),
                },
                c1: Fp2 {
                    c0: Fp(0x2b486cddf176f3ce),
                    c1: Fp(0x35208dc75e7be2c0),
                },
                c2: Fp2 {
                    c0: Fp(0x32fe2d6195af36e1),
                    c1: Fp(0x2114b20fee12d330),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x29a58de990080082),
                    c1: Fp(0x13dd498a04234b89),
                },
                c1: Fp2 {
                    c0: Fp(0x3a05eb9f5fbb07f8),
                    c1: Fp(0x2a123e4206c9690c),
                },
                c2: Fp2 {
                    c0: Fp(0x1b7e8eaadb27e62d),
                    c1: Fp(0x3381b02b76d410c0),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x36685937115ce890),
                    c1: Fp(0x2e5908916e561856),
                },
                c1: Fp2 {
                    c0: Fp(0x17c98bbd5767cdab),
                    c1: Fp(0x3ebff75754c67c0b),
                },
                c2: Fp2 {
                    c0: Fp(0xa2513695e0610db),
                    c1: Fp(0x2bb1e4ca4d411777),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1a85d0df1d63a674),
                    c1: Fp(0xf26c72270adcc40),
                },
                c1: Fp2 {
                    c0: Fp(0x9750c1bc2b52516),
                    c1: Fp(0x2167534e39881f4c),
                },
                c2: Fp2 {
                    c0: Fp(0x3690182c16c8dbf),
                    c1: Fp(0x1b3c2d37cc1079ec),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3e1158e9b045c86e),
                    c1: Fp(0x11d0be2bfe364007),
                },
                c1: Fp2 {
                    c0: Fp(0x1075ef9213b896f8),
                    c1: Fp(0x3b31bc3d5d92d3d8),
                },
                c2: Fp2 {
                    c0: Fp(0x30db07ed8ad4535f),
                    c1: Fp(0x3ff64674ca5d76d6),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xf6cf19398f54888),
                    c1: Fp(0x4dacf29cd4b5c77),
                },
                c1: Fp2 {
                    c0: Fp(0x16cb81e96ff0c648),
                    c1: Fp(0x2217eafa7d1ce4bc),
                },
                c2: Fp2 {
                    c0: Fp(0xb292440af2dfee9),
                    c1: Fp(0x34e0b09ee6f7c776),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x210a551f999e61f9),
                    c1: Fp(0x377ea7139980b630),
                },
                c1: Fp2 {
                    c0: Fp(0x3f23c89d0c4a696d),
                    c1: Fp(0x2e3155a20d8ff554),
                },
                c2: Fp2 {
                    c0: Fp(0x1be0ce017fe7c841),
                    c1: Fp(0xc29754bbb721),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x943d5b444946481),
                    c1: Fp(0xc4a8bfc03f3ad70),
                },
                c1: Fp2 {
                    c0: Fp(0x71ed5d7f898fbae),
                    c1: Fp(0x39e9183cb124651d),
                },
                c2: Fp2 {
                    c0: Fp(0x33091acf37dd9815),
                    c1: Fp(0x2f6897e77005bd2b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2e3aae8aab2f8e25),
                    c1: Fp(0x1a51646b9c60e8d1),
                },
                c1: Fp2 {
                    c0: Fp(0x16c631e714f55b5d),
                    c1: Fp(0x3ed93648563c6ed4),
                },
                c2: Fp2 {
                    c0: Fp(0x4b40489a3aef6af),
                    c1: Fp(0x3576fd368ea6d855),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x49b6f5cf7c91b7d),
                    c1: Fp(0x2a39a440c3c72733),
                },
                c1: Fp2 {
                    c0: Fp(0x15f0cf8e93c31c64),
                    c1: Fp(0x88722da3dee1583),
                },
                c2: Fp2 {
                    c0: Fp(0x2eca90952382eeb2),
                    c1: Fp(0x252da352bbfe399b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2151478922a02160),
                    c1: Fp(0x16b2af298624965a),
                },
                c1: Fp2 {
                    c0: Fp(0x3971ef9997a2513a),
                    c1: Fp(0x899be5d865f1aa4),
                },
                c2: Fp2 {
                    c0: Fp(0x17f8fc869be1f1e),
                    c1: Fp(0x233a7e5042de47d8),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ee6608520a06c2e),
                    c1: Fp(0xea7fbaea450b87d),
                },
                c1: Fp2 {
                    c0: Fp(0x3862cf3b5b20026d),
                    c1: Fp(0x1db76ae694b3a4be),
                },
                c2: Fp2 {
                    c0: Fp(0x31f83c05b67a5b67),
                    c1: Fp(0x1ca810624fe7a2ba),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ef9ccde998a6219),
                    c1: Fp(0x3d74f272ebc33e49),
                },
                c1: Fp2 {
                    c0: Fp(0x7fb9cb8721a76c1),
                    c1: Fp(0x2da3ba8c5e9e2db5),
                },
                c2: Fp2 {
                    c0: Fp(0xd0f92ccd2511936),
                    c1: Fp(0x15c9a2b81d7860c2),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25e2644d0c81e36b),
                    c1: Fp(0x30441189b7bedd15),
                },
                c1: Fp2 {
                    c0: Fp(0x20c1f9bf464db68e),
                    c1: Fp(0xd3ab36ae88f92d3),
                },
                c2: Fp2 {
                    c0: Fp(0x29893867ed390314),
                    c1: Fp(0x854714c5e8d56a6),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2a07e637cf64a4d9),
                    c1: Fp(0x19d24f70d6158c7f),
                },
                c1: Fp2 {
                    c0: Fp(0x3e12a30dafeb943e),
                    c1: Fp(0x2dfdc0991a317f8a),
                },
                c2: Fp2 {
                    c0: Fp(0x3a03f9264b9a6d3d),
                    c1: Fp(0xe4ed56d5a54dd61),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x17964f806da554d3),
                    c1: Fp(0x390d9cd3d71f54f9),
                },
                c1: Fp2 {
                    c0: Fp(0xf7523d31d125e25),
                    c1: Fp(0x14f01895b9beb379),
                },
                c2: Fp2 {
                    c0: Fp(0x1195c56084a59d16),
                    c1: Fp(0x3c6602cd8dcdddf8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x13e9590ed1785d86),
                    c1: Fp(0x317ecbcb51d21952),
                },
                c1: Fp2 {
                    c0: Fp(0x151e6d544af4cef3),
                    c1: Fp(0x2c29f07cb56923cc),
                },
                c2: Fp2 {
                    c0: Fp(0xd5053b121d87953),
                    c1: Fp(0x97e3aff3d853663),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2397e845d0ca4816),
                    c1: Fp(0x204a412409821c1),
                },
                c1: Fp2 {
                    c0: Fp(0x17e5e3b84ff0eb55),
                    c1: Fp(0x32f3c65757db41d0),
                },
                c2: Fp2 {
                    c0: Fp(0x354e3aafc2bc8c45),
                    c1: Fp(0x3828148d9b584d73),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x4169881a81f3dae4),
                    c1: Fp(0x2a27cb7b449f6567),
                },
                c1: Fp2 {
                    c0: Fp(0x12fdf96b3c69d274),
                    c1: Fp(0x2a1506a3ec9b0e7),
                },
                c2: Fp2 {
                    c0: Fp(0x16dcb5f2cf88f22e),
                    c1: Fp(0x3a64213a0ee14991),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xfb14d6c2b7d67b9),
                    c1: Fp(0x3e5ab00549220ae2),
                },
                c1: Fp2 {
                    c0: Fp(0xe53ac6a4018c399),
                    c1: Fp(0x16bca6b4fa61fd7c),
                },
                c2: Fp2 {
                    c0: Fp(0x28a09c4d611d7516),
                    c1: Fp(0x2d690a0f396d1dc3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x19467a7c68f11909),
                    c1: Fp(0xc1f5315df6b9a49),
                },
                c1: Fp2 {
                    c0: Fp(0x24d8890c9b957b6d),
                    c1: Fp(0x167c03411e41d58a),
                },
                c2: Fp2 {
                    c0: Fp(0x1472693b435d9ef1),
                    c1: Fp(0x7b42f926febe83),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1b17ad97ba25bad7),
                    c1: Fp(0x1f7ffc875c0771c5),
                },
                c1: Fp2 {
                    c0: Fp(0x7146ead1c02041),
                    c1: Fp(0x3a9cfb8a35671b96),
                },
                c2: Fp2 {
                    c0: Fp(0x1a06f09b63ef8d25),
                    c1: Fp(0x10d38cf04f746665),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xdd9f9121cae9f73),
                    c1: Fp(0x20ab9652af51b58a),
                },
                c1: Fp2 {
                    c0: Fp(0x3a96912264ab79d2),
                    c1: Fp(0x29da1e720bf13e65),
                },
                c2: Fp2 {
                    c0: Fp(0x17853c0705a57367),
                    c1: Fp(0x3686d6d13f3b213a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3dd24262227de355),
                    c1: Fp(0x1b35b45efae697ac),
                },
                c1: Fp2 {
                    c0: Fp(0x3a29794044131f7e),
                    c1: Fp(0x22c00eb34e58e13e),
                },
                c2: Fp2 {
                    c0: Fp(0x24cfedf6167c3d53),
                    c1: Fp(0x3cc1fb7ee61146a7),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1970973731cdb6e4),
                    c1: Fp(0x1c7b9667e221c4f0),
                },
                c1: Fp2 {
                    c0: Fp(0x24ab0e9343e9fef7),
                    c1: Fp(0x3cb01a87ddb0c128),
                },
                c2: Fp2 {
                    c0: Fp(0x1bf37b4ad2fa44ad),
                    c1: Fp(0x37b6bfb9df4c1b47),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xc36a3e9a7ca8c4a),
                    c1: Fp(0x30f007a2f78fb29c),
                },
                c1: Fp2 {
                    c0: Fp(0x1c0825dfab0a61b),
                    c1: Fp(0x13f6a20df6cf49ae),
                },
                c2: Fp2 {
                    c0: Fp(0x3ab422891f77f866),
                    c1: Fp(0x117597cd2633e723),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3e4f77c12f077419),
                    c1: Fp(0x300b0c60aaa01f5c),
                },
                c1: Fp2 {
                    c0: Fp(0xe5c567c9eaef2c6),
                    c1: Fp(0x3a8e9e88a8c401f9),
                },
                c2: Fp2 {
                    c0: Fp(0x26e0597004e70647),
                    c1: Fp(0x3bc89467a15af24b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2697e3250b556030),
                    c1: Fp(0x5d7fa93e03b56d0),
                },
                c1: Fp2 {
                    c0: Fp(0x165cc2ffd6f38499),
                    c1: Fp(0xdcaed4cee40c811),
                },
                c2: Fp2 {
                    c0: Fp(0x3078d6b706c39ba4),
                    c1: Fp(0x1ef0226df4d96c55),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3f462aa3e3893dac),
                    c1: Fp(0x3e5cf77696cbc05c),
                },
                c1: Fp2 {
                    c0: Fp(0x25cf5ee660d66c1a),
                    c1: Fp(0x290563293ea88c1c),
                },
                c2: Fp2 {
                    c0: Fp(0x39de54675b8fede7),
                    c1: Fp(0x34b2a28c651fe41d),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x20e2e0b571591649),
                    c1: Fp(0x3cba1ec6c1d01cad),
                },
                c1: Fp2 {
                    c0: Fp(0x3e8f4f0afd0787b),
                    c1: Fp(0x29025de34c8e9d6d),
                },
                c2: Fp2 {
                    c0: Fp(0x26a717b71780f88a),
                    c1: Fp(0x2d89d039f0a5fc6),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2af15a085bdfa7d8),
                    c1: Fp(0x5254b63e7b61573),
                },
                c1: Fp2 {
                    c0: Fp(0xd7effd4b5814a59),
                    c1: Fp(0x1ed7b2961ec573ee),
                },
                c2: Fp2 {
                    c0: Fp(0x553db4fe77505ef),
                    c1: Fp(0x1dbe85b6c918e2d7),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2fb4bd16ccaafe3d),
                    c1: Fp(0x4c6294c1477333e),
                },
                c1: Fp2 {
                    c0: Fp(0x8c2215e5db36f6f),
                    c1: Fp(0x3ec7504ee261b1dc),
                },
                c2: Fp2 {
                    c0: Fp(0xb2e3caeb26ee0cb),
                    c1: Fp(0x322928ebbf7f61d0),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ffefba7611b9de9),
                    c1: Fp(0x206073968d4bf179),
                },
                c1: Fp2 {
                    c0: Fp(0x3828e7084e13a327),
                    c1: Fp(0x24ac52a75de41a47),
                },
                c2: Fp2 {
                    c0: Fp(0x231bd1c6bdacfda7),
                    c1: Fp(0xc04d4e7ef634c64),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x17d831551279cf5a),
                    c1: Fp(0x18e45b6d3707b6f7),
                },
                c1: Fp2 {
                    c0: Fp(0x28bd25f0463ce638),
                    c1: Fp(0x3601f8bdfe83cc3a),
                },
                c2: Fp2 {
                    c0: Fp(0x24389e82df9bd917),
                    c1: Fp(0x33098d96db04f1f9),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x123fdc4758d51f94),
                    c1: Fp(0xce23291933dc2e),
                },
                c1: Fp2 {
                    c0: Fp(0x2eba3d6400585cd5),
                    c1: Fp(0x3a11efcac410766f),
                },
                c2: Fp2 {
                    c0: Fp(0x40c0ab2d93bc2e5),
                    c1: Fp(0x6ca9550485312ad),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x266dea4e1f4c99c4),
                    c1: Fp(0x30820d6d89e0b128),
                },
                c1: Fp2 {
                    c0: Fp(0x30f91076061f10f7),
                    c1: Fp(0x25032007b70ae2f7),
                },
                c2: Fp2 {
                    c0: Fp(0x2b92cc82eda0fd44),
                    c1: Fp(0x11c2e96ae55850d4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1151b3e025c03d10),
                    c1: Fp(0x2ae6b0fd6c40841d),
                },
                c1: Fp2 {
                    c0: Fp(0x28fed8b79338bf07),
                    c1: Fp(0x2a0fdbab19cff68),
                },
                c2: Fp2 {
                    c0: Fp(0x394844db5ead9291),
                    c1: Fp(0x15c8fbbcba143022),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28b3b4df996afb2b),
                    c1: Fp(0x1449a512c0e5b812),
                },
                c1: Fp2 {
                    c0: Fp(0xec5ea61bb4ff354),
                    c1: Fp(0x26e20e7609861980),
                },
                c2: Fp2 {
                    c0: Fp(0x2beef1a0be865c50),
                    c1: Fp(0x1dccd3334e08dd94),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1d09fc7a0f20951e),
                    c1: Fp(0x1f2fc573acdeec25),
                },
                c1: Fp2 {
                    c0: Fp(0x3d5d2bad79fa37fe),
                    c1: Fp(0x292da96c1f8ec493),
                },
                c2: Fp2 {
                    c0: Fp(0x161cdd1157a88197),
                    c1: Fp(0x2fc8e11908e7433f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xbaa6397a02b7a9c),
                    c1: Fp(0x1edb16e74ed05c4e),
                },
                c1: Fp2 {
                    c0: Fp(0x2ca486567931c395),
                    c1: Fp(0x3831cce9e17f22f7),
                },
                c2: Fp2 {
                    c0: Fp(0x128761df4a1de088),
                    c1: Fp(0x4118ee5984e59be2),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x35106bb5c584c60),
                    c1: Fp(0x4cd6d6ebdb476f3),
                },
                c1: Fp2 {
                    c0: Fp(0x2fd344d6e9dce669),
                    c1: Fp(0x40cfb9a6ebb80948),
                },
                c2: Fp2 {
                    c0: Fp(0x3a200ed610a16914),
                    c1: Fp(0x30fb19f0e42e2ae3),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x279ece3c731fa4c2),
                    c1: Fp(0x3119352810b22b2b),
                },
                c1: Fp2 {
                    c0: Fp(0xc4302d84bbfc06e),
                    c1: Fp(0x24fe2eda7c2bcca8),
                },
                c2: Fp2 {
                    c0: Fp(0x33c4f14ae39994cd),
                    c1: Fp(0x1e7f108eeeafd4fd),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1d2e94793e2c735d),
                    c1: Fp(0x376c4590b8381818),
                },
                c1: Fp2 {
                    c0: Fp(0x282db3c9ac054c99),
                    c1: Fp(0x2801bb8621e5e71d),
                },
                c2: Fp2 {
                    c0: Fp(0x36e42543564c1b30),
                    c1: Fp(0x1e1bcaf06f8228d0),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x939f334723a5de2),
                    c1: Fp(0x2fcdd89d042c7d9d),
                },
                c1: Fp2 {
                    c0: Fp(0x22359c770a11f903),
                    c1: Fp(0x275b5f46d0821f66),
                },
                c2: Fp2 {
                    c0: Fp(0x408698a828cb85a3),
                    c1: Fp(0x3e1e59f170d8c857),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28245e10e53f4811),
                    c1: Fp(0x2baa604adbc56411),
                },
                c1: Fp2 {
                    c0: Fp(0x3aae92b6b46f1aa7),
                    c1: Fp(0x34ced76fdb676324),
                },
                c2: Fp2 {
                    c0: Fp(0x3380be2a45d353ae),
                    c1: Fp(0x2676ec3d22fec84e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1c47dbfd43940f5e),
                    c1: Fp(0x195070c41342465),
                },
                c1: Fp2 {
                    c0: Fp(0x1080a1d4b1780b78),
                    c1: Fp(0x8770676e382eb56),
                },
                c2: Fp2 {
                    c0: Fp(0x1d4566d52c8385b2),
                    c1: Fp(0x15df59b8e05f7c67),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a6fe7a787970f64),
                    c1: Fp(0x12e488e86509d30c),
                },
                c1: Fp2 {
                    c0: Fp(0x16374150c71e32e0),
                    c1: Fp(0x15a608633fb69a4d),
                },
                c2: Fp2 {
                    c0: Fp(0x1a0d46d2ee8eba58),
                    c1: Fp(0x1ed031933267d9a0),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x4ff5fb2e8c30ff3),
                    c1: Fp(0x340b90a2bdd869a6),
                },
                c1: Fp2 {
                    c0: Fp(0x35533ddbace957ce),
                    c1: Fp(0x30df9d9cfa62fee1),
                },
                c2: Fp2 {
                    c0: Fp(0x3a353f6bc9ec59f),
                    c1: Fp(0xca698cad8b7e39f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3eb5cb48341ce60e),
                    c1: Fp(0x3532aded18a01c2d),
                },
                c1: Fp2 {
                    c0: Fp(0x8726da93d452dac),
                    c1: Fp(0x4d11da258a402bf),
                },
                c2: Fp2 {
                    c0: Fp(0x28d84fc09a794b39),
                    c1: Fp(0x4de43acba58ecab),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd737e0e7b80c77a),
                    c1: Fp(0x20fc3bb674b1095d),
                },
                c1: Fp2 {
                    c0: Fp(0x3229f62e8ee6e37),
                    c1: Fp(0x372502594877bc5c),
                },
                c2: Fp2 {
                    c0: Fp(0x1c368987996d2b5f),
                    c1: Fp(0x3e5416d1cf3da56f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x15737bd1556fe9ed),
                    c1: Fp(0x1fecba0616789851),
                },
                c1: Fp2 {
                    c0: Fp(0x27fbd635187a5f3),
                    c1: Fp(0xbd9bd7ca056af97),
                },
                c2: Fp2 {
                    c0: Fp(0x2da3818670b69418),
                    c1: Fp(0x1d1bbe495ed77241),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x320fc92add444717),
                    c1: Fp(0x17bbe4919230c5f),
                },
                c1: Fp2 {
                    c0: Fp(0x3e78adab3e9e8e51),
                    c1: Fp(0x37552b422799c032),
                },
                c2: Fp2 {
                    c0: Fp(0x2152b3c80b02f2d4),
                    c1: Fp(0x701b1144b2e3b8a),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3582695953925342),
                    c1: Fp(0x3672a3cdde13ad67),
                },
                c1: Fp2 {
                    c0: Fp(0x78acd050c0b32ef),
                    c1: Fp(0x1f04d5da882c0f58),
                },
                c2: Fp2 {
                    c0: Fp(0x150862a8b88725c6),
                    c1: Fp(0x8d72fba85e8e9),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3672727eea2e3f4c),
                    c1: Fp(0x2fe6e4dbc6726b56),
                },
                c1: Fp2 {
                    c0: Fp(0x296a78cc1b04c9ea),
                    c1: Fp(0x2ffb00eeb6d04a71),
                },
                c2: Fp2 {
                    c0: Fp(0x8a1b98e94637620),
                    c1: Fp(0xc426d4fb81fe62d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3b5cdb2c20de631d),
                    c1: Fp(0x40c4871c8d0c8762),
                },
                c1: Fp2 {
                    c0: Fp(0x2e7243601c8d0e97),
                    c1: Fp(0x2c7981b29f7c1cf6),
                },
                c2: Fp2 {
                    c0: Fp(0x36aa53701b7703f6),
                    c1: Fp(0x344e017ea23436f5),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xa59c4d74259f933),
                    c1: Fp(0x3800016a23cd9d11),
                },
                c1: Fp2 {
                    c0: Fp(0x1cca12b835bff62e),
                    c1: Fp(0x2261d2b88f19edf2),
                },
                c2: Fp2 {
                    c0: Fp(0x1c5e9163c8676e7a),
                    c1: Fp(0x2717ba06be68fa84),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x169850fa19685743),
                    c1: Fp(0x16b7a502ccc74d89),
                },
                c1: Fp2 {
                    c0: Fp(0x1faa8eac48f31ea7),
                    c1: Fp(0x1c8331ff5c800b1d),
                },
                c2: Fp2 {
                    c0: Fp(0x1a8616d269409b39),
                    c1: Fp(0x1e4fe1832a0b427b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1091c8dc968b8ff),
                    c1: Fp(0x2ac3212dc9a71ce2),
                },
                c1: Fp2 {
                    c0: Fp(0x2716bfb193d80236),
                    c1: Fp(0x1e49bc101c753b01),
                },
                c2: Fp2 {
                    c0: Fp(0x255213a61cee582f),
                    c1: Fp(0x18a7eed73bf9dcc8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ff5b8faea87a3d7),
                    c1: Fp(0x2813aa58c1f052e7),
                },
                c1: Fp2 {
                    c0: Fp(0x3a76eff9068d0e47),
                    c1: Fp(0x2a10f408347d8aa4),
                },
                c2: Fp2 {
                    c0: Fp(0xdc354dc73c87a61),
                    c1: Fp(0x25c52b0c931f062e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x37c8db35cea25d95),
                    c1: Fp(0x3d2045cd327cde22),
                },
                c1: Fp2 {
                    c0: Fp(0x2a60ead1d2d8b521),
                    c1: Fp(0x2826889c2d6989e1),
                },
                c2: Fp2 {
                    c0: Fp(0x3b96c7a93d2e22b6),
                    c1: Fp(0x6212c742c23150c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25f0d5220ca7f328),
                    c1: Fp(0x1afb7d72e5d3d1ef),
                },
                c1: Fp2 {
                    c0: Fp(0x182003cc1bb7d1e),
                    c1: Fp(0x486a57355077f03),
                },
                c2: Fp2 {
                    c0: Fp(0xe035bc3d6c80104),
                    c1: Fp(0x1cb4ff867134b8be),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2d1accc11e88778a),
                    c1: Fp(0x5503208f5a375a1),
                },
                c1: Fp2 {
                    c0: Fp(0x2da7fe7703c43916),
                    c1: Fp(0x155dd6225d462db5),
                },
                c2: Fp2 {
                    c0: Fp(0x38b4fef8554fb78c),
                    c1: Fp(0x30e716f95eb1e3b8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x23ec9c2ffa2de145),
                    c1: Fp(0x256712e99c69720c),
                },
                c1: Fp2 {
                    c0: Fp(0x2c1d97203cb44780),
                    c1: Fp(0x2e92040de9d1ce7e),
                },
                c2: Fp2 {
                    c0: Fp(0x233d67af1d3b5646),
                    c1: Fp(0x16e0f94e4dbba0f7),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x235b2e1a1e50a474),
                    c1: Fp(0x1dcb07d7e4f0be4d),
                },
                c1: Fp2 {
                    c0: Fp(0x1d8891b851627118),
                    c1: Fp(0xb3496fd5c8d85e2),
                },
                c2: Fp2 {
                    c0: Fp(0x3c12ca6b07b42635),
                    c1: Fp(0x64e792fa7fa82b1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1b6459992ea54c28),
                    c1: Fp(0x29d9bb457ed579c8),
                },
                c1: Fp2 {
                    c0: Fp(0x2bb6d345e3affe41),
                    c1: Fp(0x29ba055fe1961881),
                },
                c2: Fp2 {
                    c0: Fp(0x3591a9f10ae2c8d9),
                    c1: Fp(0xfa3bd633d003296),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x4118f87e47fad35d),
                    c1: Fp(0x37845ce834769ef2),
                },
                c1: Fp2 {
                    c0: Fp(0x2b786603f3ef066c),
                    c1: Fp(0x19838c60991931),
                },
                c2: Fp2 {
                    c0: Fp(0x1aa8c21fb546bf9b),
                    c1: Fp(0x3ce89e0d7bba35dd),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x33bbd6c813043fcf),
                    c1: Fp(0x3b9a79170a25e484),
                },
                c1: Fp2 {
                    c0: Fp(0x2ccbe79625cec859),
                    c1: Fp(0x1c4d721674076b83),
                },
                c2: Fp2 {
                    c0: Fp(0x313b99fb00f30390),
                    c1: Fp(0x2433cbb808f58b63),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe9b6b12d3763a2c),
                    c1: Fp(0x144c27932f9acb14),
                },
                c1: Fp2 {
                    c0: Fp(0x5f65941e67172f7),
                    c1: Fp(0x3eefb1a5af6b1a02),
                },
                c2: Fp2 {
                    c0: Fp(0x40ae6b52180c1bc0),
                    c1: Fp(0x3986f7feaeebe658),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x34437b278abd2c66),
                    c1: Fp(0x2bd7b37e9cefd04c),
                },
                c1: Fp2 {
                    c0: Fp(0x368419828fe1e3d0),
                    c1: Fp(0x3eee18da0a35116f),
                },
                c2: Fp2 {
                    c0: Fp(0x244cfdb916349612),
                    c1: Fp(0x39da08a26efc70e8),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x17d327c07ab31817),
                    c1: Fp(0x9e5006772ec2908),
                },
                c1: Fp2 {
                    c0: Fp(0x63b2b5196f1db98),
                    c1: Fp(0x16f715b68eab63fb),
                },
                c2: Fp2 {
                    c0: Fp(0x31f6ea981d9bcdf9),
                    c1: Fp(0x65dcbf96ffc1af4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x286bcaaa56380c78),
                    c1: Fp(0xdfc6cb6553f6470),
                },
                c1: Fp2 {
                    c0: Fp(0x3460bec1ee5dae48),
                    c1: Fp(0x20b04adfbb44b91c),
                },
                c2: Fp2 {
                    c0: Fp(0x3b7b177ddc5c492f),
                    c1: Fp(0xb8bf49cb4366195),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xa7731d325abe684),
                    c1: Fp(0x1193081fb7a9e7c1),
                },
                c1: Fp2 {
                    c0: Fp(0x3a4b7d2385f6d725),
                    c1: Fp(0xb98399010bb17ef),
                },
                c2: Fp2 {
                    c0: Fp(0x16d3550087d84944),
                    c1: Fp(0xc95eda685b60619),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3f2977ba548a4eb),
                    c1: Fp(0x1780a68bf5ec750e),
                },
                c1: Fp2 {
                    c0: Fp(0x28d1e3630166f782),
                    c1: Fp(0x277f877eb7241f8a),
                },
                c2: Fp2 {
                    c0: Fp(0x1a42e41469d62139),
                    c1: Fp(0x394765e8e8f8a19f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1a6df351fb5bffb3),
                    c1: Fp(0xad3c3b4076076c3),
                },
                c1: Fp2 {
                    c0: Fp(0x3dc7bf6f0fa6a027),
                    c1: Fp(0x2afe8326f4ad9374),
                },
                c2: Fp2 {
                    c0: Fp(0x14b2c3b7c6398136),
                    c1: Fp(0xfa53732d509a095),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x31c43cdc3537c82c),
                    c1: Fp(0xd8b87b2a8ee47ae),
                },
                c1: Fp2 {
                    c0: Fp(0x35d20886e2937ce1),
                    c1: Fp(0x20f5fcbd89dfdff2),
                },
                c2: Fp2 {
                    c0: Fp(0x145a4d462db77756),
                    c1: Fp(0x25cafb3bb198aa89),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x674ecb7f6aafcfa),
                    c1: Fp(0x3015796da9b99906),
                },
                c1: Fp2 {
                    c0: Fp(0x12b6dd54475ada80),
                    c1: Fp(0x32b607b1193e03f),
                },
                c2: Fp2 {
                    c0: Fp(0x3cd091dceb48b903),
                    c1: Fp(0x6fe9e0f8caeb878),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3bd0a2b7547846a),
                    c1: Fp(0x3dba673e5c3631d5),
                },
                c1: Fp2 {
                    c0: Fp(0x230eca9769e10e9a),
                    c1: Fp(0x234ccd0f6a81859a),
                },
                c2: Fp2 {
                    c0: Fp(0x3ef99a2ebb7cf891),
                    c1: Fp(0x3b5cf2de6cc9bcea),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1266f0024125b400),
                    c1: Fp(0xc2f35acb7844d75),
                },
                c1: Fp2 {
                    c0: Fp(0x16e1668a7336512e),
                    c1: Fp(0x2781caef8eebe120),
                },
                c2: Fp2 {
                    c0: Fp(0x3706e26b8e80032d),
                    c1: Fp(0x20475391ad5f235e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd2ac4218f3da38b),
                    c1: Fp(0x1e93e79ef34829a1),
                },
                c1: Fp2 {
                    c0: Fp(0x3b13f60092b0e9b2),
                    c1: Fp(0x15ceb34f26e00815),
                },
                c2: Fp2 {
                    c0: Fp(0xac684e3ac0e48ab),
                    c1: Fp(0x1e27ace28ac60fca),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x30213649fa7cee16),
                    c1: Fp(0x2696b4216577d14b),
                },
                c1: Fp2 {
                    c0: Fp(0x17cf87cd9266e619),
                    c1: Fp(0x10d9821ffa0f62f6),
                },
                c2: Fp2 {
                    c0: Fp(0x1559e7cbd39b3ecd),
                    c1: Fp(0x2c8fe6570e3877a2),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd1fd7c6ff714f08),
                    c1: Fp(0x14dcafec109b1a93),
                },
                c1: Fp2 {
                    c0: Fp(0x1971fddf26d0df31),
                    c1: Fp(0x402e1050d23203b6),
                },
                c2: Fp2 {
                    c0: Fp(0x1d5041bc9e493ce2),
                    c1: Fp(0x2f134f05355272da),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1b96d0165cf00e69),
                    c1: Fp(0x23c5a6de8e153553),
                },
                c1: Fp2 {
                    c0: Fp(0x2783cac51e2322f4),
                    c1: Fp(0x39c1eb85e113f522),
                },
                c2: Fp2 {
                    c0: Fp(0x28530f6903a5b00a),
                    c1: Fp(0x1cfbf28ed81bb387),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c6d5f9a0e247861),
                    c1: Fp(0x25e6defcd8a53bba),
                },
                c1: Fp2 {
                    c0: Fp(0x6fe4fb445c78e84),
                    c1: Fp(0x343c36af674da4e2),
                },
                c2: Fp2 {
                    c0: Fp(0x3951e424d04e245b),
                    c1: Fp(0x12b3aca7791d1ed2),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x29efc6d5b9c0f0a8),
                    c1: Fp(0x1e81d2d86001f0c),
                },
                c1: Fp2 {
                    c0: Fp(0x3482ccb9d6525897),
                    c1: Fp(0x10544fd5c8c7e58b),
                },
                c2: Fp2 {
                    c0: Fp(0x3d7ab4552caf7169),
                    c1: Fp(0x3ac288d12eebf49c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x40c09d8c1eb4d444),
                    c1: Fp(0x18a6b7f8a49184b8),
                },
                c1: Fp2 {
                    c0: Fp(0x2ec430cd5f5639d2),
                    c1: Fp(0x2bdc44e9a799f5d0),
                },
                c2: Fp2 {
                    c0: Fp(0x2e2e1e0247fd4cee),
                    c1: Fp(0xa5db6d56a7d83e3),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x48712e2cd4d5688),
                    c1: Fp(0x2dd4a01b3f7f9557),
                },
                c1: Fp2 {
                    c0: Fp(0x117db0b1e16d65a9),
                    c1: Fp(0x144462f08b661b61),
                },
                c2: Fp2 {
                    c0: Fp(0xd7827b1df4bf60d),
                    c1: Fp(0x1bbbdf3485c6e93e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x35f90d88e190d17),
                    c1: Fp(0x25f83cc756301ea6),
                },
                c1: Fp2 {
                    c0: Fp(0x2bd6da0430e90832),
                    c1: Fp(0x6aa65703273bb2a),
                },
                c2: Fp2 {
                    c0: Fp(0x2a72e4f05785e6),
                    c1: Fp(0x1658e3d44859628b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x8f2d1d6bdfd8eca),
                    c1: Fp(0x5e97618bb73c167),
                },
                c1: Fp2 {
                    c0: Fp(0x22232a3ea92f6543),
                    c1: Fp(0xfc1d4550fc6ede9),
                },
                c2: Fp2 {
                    c0: Fp(0x29d338c53b0d32e4),
                    c1: Fp(0x9c3f6334972a7ed),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x13714bafd118c4f6),
                    c1: Fp(0x151404481a6706ec),
                },
                c1: Fp2 {
                    c0: Fp(0x1da2cc09988643b7),
                    c1: Fp(0x3424f25aa80b1d01),
                },
                c2: Fp2 {
                    c0: Fp(0x38b31c50531391fc),
                    c1: Fp(0x12571c3d79b23814),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x412cea24f55b1b90),
                    c1: Fp(0x38e0e0998ec23e7c),
                },
                c1: Fp2 {
                    c0: Fp(0x29fe1cfb8877b160),
                    c1: Fp(0x1e5a45f048cb504e),
                },
                c2: Fp2 {
                    c0: Fp(0x33cd9397d5e2bcf1),
                    c1: Fp(0x3baca0be86ce6bde),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1fcec5c26059ce10),
                    c1: Fp(0x317819986be7047),
                },
                c1: Fp2 {
                    c0: Fp(0x3a11ee896a6d16c9),
                    c1: Fp(0x3af9b96fb4a1ead8),
                },
                c2: Fp2 {
                    c0: Fp(0x2e2494eb7ac6dcd9),
                    c1: Fp(0x1646846eabfae155),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3f022e8949259ebb),
                    c1: Fp(0x3c1257af2cf47695),
                },
                c1: Fp2 {
                    c0: Fp(0x3bde598eb60994fa),
                    c1: Fp(0x3b3682121cab309),
                },
                c2: Fp2 {
                    c0: Fp(0x1ada5463a22cdc2a),
                    c1: Fp(0xbdb1cdfbab054a8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x17f7ae0d2899c057),
                    c1: Fp(0x312902708dfa187),
                },
                c1: Fp2 {
                    c0: Fp(0x2e8b7877f1289de9),
                    c1: Fp(0x2e18563f57019b1a),
                },
                c2: Fp2 {
                    c0: Fp(0x3e3b1f0e27c22536),
                    c1: Fp(0x2811058ed39989da),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x379a72ae64f630f3),
                    c1: Fp(0x256ca157e274e717),
                },
                c1: Fp2 {
                    c0: Fp(0xe1d5df21642ad2e),
                    c1: Fp(0x2d0b04e15b8dcb8),
                },
                c2: Fp2 {
                    c0: Fp(0x1d7040a6e5d26c96),
                    c1: Fp(0x38d9026ae5cb799b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28649616d6f072fb),
                    c1: Fp(0x2c41b76238c0f967),
                },
                c1: Fp2 {
                    c0: Fp(0x73d7682032d154a),
                    c1: Fp(0x20aa64ff109a88ef),
                },
                c2: Fp2 {
                    c0: Fp(0x5eb461c616d7c39),
                    c1: Fp(0x39f55ba5cc6f1bd3),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3199040eab64f700),
                    c1: Fp(0x1c67705995d9247d),
                },
                c1: Fp2 {
                    c0: Fp(0x3aa56a110ccd3d0a),
                    c1: Fp(0x2d05a49af05ad474),
                },
                c2: Fp2 {
                    c0: Fp(0x5c76b5024113609),
                    c1: Fp(0xa31896ee3b1e836),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3571e660b1c3252e),
                    c1: Fp(0x357af3fd94b8df6c),
                },
                c1: Fp2 {
                    c0: Fp(0x2b11b5374b7fbc93),
                    c1: Fp(0x11fa5cb4d8a5c710),
                },
                c2: Fp2 {
                    c0: Fp(0x36399a755133e553),
                    c1: Fp(0x2c503910e6ee504a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3410ac2429b3389),
                    c1: Fp(0x27c91aa80d48e721),
                },
                c1: Fp2 {
                    c0: Fp(0x3d178ab4fc64380d),
                    c1: Fp(0x20328cc079819eef),
                },
                c2: Fp2 {
                    c0: Fp(0x35c7ed4abff10af0),
                    c1: Fp(0x32db1b34c6e5bce4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x37c488ded3e1ca49),
                    c1: Fp(0x19a363a08837bf7c),
                },
                c1: Fp2 {
                    c0: Fp(0x334c1e638ed46e7),
                    c1: Fp(0x3a1cd480fe2582f0),
                },
                c2: Fp2 {
                    c0: Fp(0x2883ade49461f016),
                    c1: Fp(0x2f5beafd6caf5b3c),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xcbad61a3747f1c5),
                    c1: Fp(0xe92783a676f0576),
                },
                c1: Fp2 {
                    c0: Fp(0x4086a138c95044bd),
                    c1: Fp(0xac02bc45c5bdd1),
                },
                c2: Fp2 {
                    c0: Fp(0x375f7efe914c115a),
                    c1: Fp(0x17e72a1d9a22be62),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1d567c20da6324c3),
                    c1: Fp(0x5571bfaed1ca48a),
                },
                c1: Fp2 {
                    c0: Fp(0x415888ffc5b8591a),
                    c1: Fp(0x3fd588d564fd9a5b),
                },
                c2: Fp2 {
                    c0: Fp(0x3afea2ac950563d7),
                    c1: Fp(0x1559e9f806820298),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x47b761b526591a),
                    c1: Fp(0x34d071a2071006c0),
                },
                c1: Fp2 {
                    c0: Fp(0x292d05e5de049ef6),
                    c1: Fp(0x35a7885650978e89),
                },
                c2: Fp2 {
                    c0: Fp(0x241162f955ed73e2),
                    c1: Fp(0x27b7ac69f4cb1b4f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2c4fcfaf2b31d075),
                    c1: Fp(0x3bc0165a172aa917),
                },
                c1: Fp2 {
                    c0: Fp(0x3c3702acfd0fc921),
                    c1: Fp(0x2b47a6d9ed3e0b97),
                },
                c2: Fp2 {
                    c0: Fp(0x22c84d61d2644759),
                    c1: Fp(0x5c28a479b2ff91b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1b3acd0400865aed),
                    c1: Fp(0xe47253ca19d861e),
                },
                c1: Fp2 {
                    c0: Fp(0x2a39e8c1a5c7765b),
                    c1: Fp(0x28560d1aa00722b2),
                },
                c2: Fp2 {
                    c0: Fp(0x264193c9e75239da),
                    c1: Fp(0x1a5088158743e445),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x892d56a02393b3a),
                    c1: Fp(0x22ac24db36ee18de),
                },
                c1: Fp2 {
                    c0: Fp(0x1267063834e3a9d1),
                    c1: Fp(0x22ff952adf36d2bc),
                },
                c2: Fp2 {
                    c0: Fp(0x3c11f4fd974b92d1),
                    c1: Fp(0x3e94b1e61150de02),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c246766d7271127),
                    c1: Fp(0x1ea6f2c7b1ba7eaa),
                },
                c1: Fp2 {
                    c0: Fp(0x122ca2b8e6e2e33c),
                    c1: Fp(0x351d5bfb5074134c),
                },
                c2: Fp2 {
                    c0: Fp(0x113843ac958995b3),
                    c1: Fp(0x4a0a013cbdbf81c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x5cd2f99b33f0381),
                    c1: Fp(0x3d75d854fbf133ab),
                },
                c1: Fp2 {
                    c0: Fp(0x29b5b0f23bad1474),
                    c1: Fp(0x7e2c16756be6cf2),
                },
                c2: Fp2 {
                    c0: Fp(0x305009ef75fec960),
                    c1: Fp(0x14d4da04bf855ba4),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe10eac97cb16922),
                    c1: Fp(0x244f3254b1b81419),
                },
                c1: Fp2 {
                    c0: Fp(0x4005f5dcc39ec917),
                    c1: Fp(0x116116d7baf1cba6),
                },
                c2: Fp2 {
                    c0: Fp(0x12db56aab48e7183),
                    c1: Fp(0x38b6d6dca05efb7c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x193c9a72190a380a),
                    c1: Fp(0x309958838a4bf926),
                },
                c1: Fp2 {
                    c0: Fp(0x2696544ddfc21e2b),
                    c1: Fp(0xda5931905e0ef9d),
                },
                c2: Fp2 {
                    c0: Fp(0x158b6afc26a5cdd0),
                    c1: Fp(0x28f7b8f575f03f24),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x138c64e067c84710),
                    c1: Fp(0x244f80a033ee859b),
                },
                c1: Fp2 {
                    c0: Fp(0x1a47a5143700c7c5),
                    c1: Fp(0xe2e3385bd024013),
                },
                c2: Fp2 {
                    c0: Fp(0x164b5a44cf388270),
                    c1: Fp(0x318d0aed18b9f8be),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x9e556a72bc76369),
                    c1: Fp(0x1d19c1ad5dc9a284),
                },
                c1: Fp2 {
                    c0: Fp(0x2a022254fbd3aca1),
                    c1: Fp(0x3c3153540ef43c06),
                },
                c2: Fp2 {
                    c0: Fp(0x33e2fc2ce7a71d89),
                    c1: Fp(0x1a8ff91fe3ab9c7c),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x184d234e398a929e),
                    c1: Fp(0x4e7115428ccedd5),
                },
                c1: Fp2 {
                    c0: Fp(0x1b742cd4cebc2187),
                    c1: Fp(0x2939250a10b29893),
                },
                c2: Fp2 {
                    c0: Fp(0x2c65a3d31c912152),
                    c1: Fp(0x350afa48a0e19208),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x183e4cd41fa0ab67),
                    c1: Fp(0x82b530b8ac15ac),
                },
                c1: Fp2 {
                    c0: Fp(0x1966937b52aa6b95),
                    c1: Fp(0x113a87848fcdfe1),
                },
                c2: Fp2 {
                    c0: Fp(0x23c2a25704b19afb),
                    c1: Fp(0x19a4f7b0009bd399),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x17212e39549d0662),
                    c1: Fp(0x36e2335722230e9b),
                },
                c1: Fp2 {
                    c0: Fp(0x33e399cc52bf049),
                    c1: Fp(0x153bafc9d53e7e9c),
                },
                c2: Fp2 {
                    c0: Fp(0x3c913754056595),
                    c1: Fp(0x37ee9da99a0f0fbf),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x23f93d5990ccc38f),
                    c1: Fp(0x39315f80594b895a),
                },
                c1: Fp2 {
                    c0: Fp(0x411159f03f84deb4),
                    c1: Fp(0x120936d1b229bf46),
                },
                c2: Fp2 {
                    c0: Fp(0x2f92335eb7752880),
                    c1: Fp(0x262258272156bbfa),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1fe604e53b1618bc),
                    c1: Fp(0x15c25d4f69454a66),
                },
                c1: Fp2 {
                    c0: Fp(0x345519e1ea3fb422),
                    c1: Fp(0x3878111ed391eb72),
                },
                c2: Fp2 {
                    c0: Fp(0x142eaade50a3b6d4),
                    c1: Fp(0x32e10970347596aa),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2accc4424b7690ea),
                    c1: Fp(0x10de80d7049a7964),
                },
                c1: Fp2 {
                    c0: Fp(0x20654776777bfecb),
                    c1: Fp(0x1aa712d7d95860e5),
                },
                c2: Fp2 {
                    c0: Fp(0x28e32be25bdd62c2),
                    c1: Fp(0x40c4f085a7d99c13),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x31f8d0428e5bbc8),
                    c1: Fp(0xd26d627a47c1ebf),
                },
                c1: Fp2 {
                    c0: Fp(0x8b2b0621bb03030),
                    c1: Fp(0x1f31f31b57e942cd),
                },
                c2: Fp2 {
                    c0: Fp(0x29caed7379a8231a),
                    c1: Fp(0x196b141b3eef2c6f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xfb176c189f3b510),
                    c1: Fp(0x3d86b3f0d180ba03),
                },
                c1: Fp2 {
                    c0: Fp(0x3ce5543993af1969),
                    c1: Fp(0xe0ce40890925964),
                },
                c2: Fp2 {
                    c0: Fp(0x83a79e72cf782dd),
                    c1: Fp(0xf067c0b161bd458),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x26b60723a0e2416a),
                    c1: Fp(0x23f8d72f69ccfc5c),
                },
                c1: Fp2 {
                    c0: Fp(0x1aab273331f79121),
                    c1: Fp(0x29fd41720e0b770),
                },
                c2: Fp2 {
                    c0: Fp(0x15953b8f441f60ac),
                    c1: Fp(0x2444c6c4f4bb8c71),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1017f4722c9397ed),
                    c1: Fp(0x3feefb24f8dad7f6),
                },
                c1: Fp2 {
                    c0: Fp(0x15ddde4d1ec41db),
                    c1: Fp(0x6d8f656682778a),
                },
                c2: Fp2 {
                    c0: Fp(0x3b75aa7a4f27b99f),
                    c1: Fp(0x28c48b4c2260ce1),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xefe308ffe9ed1f0),
                    c1: Fp(0x41090e100320c2e7),
                },
                c1: Fp2 {
                    c0: Fp(0x1f98defe7a37cf3a),
                    c1: Fp(0x2814b8f49d9618d7),
                },
                c2: Fp2 {
                    c0: Fp(0x154a99a5808056e9),
                    c1: Fp(0x3e816c14ebe54a8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xf4d15b09cd2007f),
                    c1: Fp(0x10c15085c3c6c60c),
                },
                c1: Fp2 {
                    c0: Fp(0x1e9b455d099d538f),
                    c1: Fp(0x1cad8448bac540d0),
                },
                c2: Fp2 {
                    c0: Fp(0x22593c0b391bbdca),
                    c1: Fp(0x17a2f9bc6aeb8acf),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2e56a899e9c574dc),
                    c1: Fp(0x200e58760f65bdd2),
                },
                c1: Fp2 {
                    c0: Fp(0x2e0df8ea348e2a9e),
                    c1: Fp(0x3ad634c67ac15873),
                },
                c2: Fp2 {
                    c0: Fp(0x173928f377ff0640),
                    c1: Fp(0xb5260c3d0c1b807),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x33afe2d0202a9899),
                    c1: Fp(0x255a5bcc08350560),
                },
                c1: Fp2 {
                    c0: Fp(0x8b7db91db33bc52),
                    c1: Fp(0x2cc395e0e4efeb77),
                },
                c2: Fp2 {
                    c0: Fp(0x2ebc64dbca103da8),
                    c1: Fp(0x3db9cdc55dd1b012),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xf146ef542729cf4),
                    c1: Fp(0x104cda0726d89f51),
                },
                c1: Fp2 {
                    c0: Fp(0x13b263bb00eb08c2),
                    c1: Fp(0x3727c23766f5ae29),
                },
                c2: Fp2 {
                    c0: Fp(0x1147c57300b7d71a),
                    c1: Fp(0x14841ac15c35a0ea),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x475ab387825ae26),
                    c1: Fp(0x40f4263660073e82),
                },
                c1: Fp2 {
                    c0: Fp(0x201e3ceefa7235e3),
                    c1: Fp(0x209211ccedea02c7),
                },
                c2: Fp2 {
                    c0: Fp(0xdf1ee4c601507d3),
                    c1: Fp(0x10babc8bce28663a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x17e28e9246d2a3bb),
                    c1: Fp(0x3c82a79eeae7238e),
                },
                c1: Fp2 {
                    c0: Fp(0x30a0d5c71fc133ef),
                    c1: Fp(0xd55e82db1ff55a),
                },
                c2: Fp2 {
                    c0: Fp(0xfcaee69dceb8ea0),
                    c1: Fp(0x3c3bf2b61eb6a7e7),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x383d8b0beb0f7764),
                    c1: Fp(0x32615a4e1a227fb1),
                },
                c1: Fp2 {
                    c0: Fp(0x40c1b7d5ad253d3b),
                    c1: Fp(0x3452c0e87a920029),
                },
                c2: Fp2 {
                    c0: Fp(0x4142ec91100892c3),
                    c1: Fp(0xf82de63212a9349),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x297d243f72013bd6),
                    c1: Fp(0x15f388a6d9d41ca0),
                },
                c1: Fp2 {
                    c0: Fp(0x2675b2051a600774),
                    c1: Fp(0x22cf185c5cad745e),
                },
                c2: Fp2 {
                    c0: Fp(0x14572766fd26ea16),
                    c1: Fp(0x2ec2039d97230762),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ac8e307fd5029c7),
                    c1: Fp(0x367ba9f27280c4fd),
                },
                c1: Fp2 {
                    c0: Fp(0x326e9d0df904180a),
                    c1: Fp(0x6af129d513f683b),
                },
                c2: Fp2 {
                    c0: Fp(0x145eba537847b3ba),
                    c1: Fp(0xa6947c1cc4de8b4),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x35ace11388d80d9f),
                    c1: Fp(0x27b6e2dc46b363b3),
                },
                c1: Fp2 {
                    c0: Fp(0x3348a485e0ad904a),
                    c1: Fp(0x268b9d098ee8d348),
                },
                c2: Fp2 {
                    c0: Fp(0x26d7b26d54936578),
                    c1: Fp(0x290a4c60e9fb3c36),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3fa7ab576ecf652d),
                    c1: Fp(0xb94f8d1180412ba),
                },
                c1: Fp2 {
                    c0: Fp(0x3cafbc36cea92220),
                    c1: Fp(0x1f47a18de37b8926),
                },
                c2: Fp2 {
                    c0: Fp(0x29e5d28d925b11fe),
                    c1: Fp(0x17da48b01eac592b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3be699cae7341033),
                    c1: Fp(0x34ba0b99d5a6439d),
                },
                c1: Fp2 {
                    c0: Fp(0xc6adfd9b7874e18),
                    c1: Fp(0x4174b5849eeaed),
                },
                c2: Fp2 {
                    c0: Fp(0x135dfa87f113363d),
                    c1: Fp(0x3b4c2f45cbbc96e3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2ee6f757fc0671f5),
                    c1: Fp(0x5c06a328bbf2046),
                },
                c1: Fp2 {
                    c0: Fp(0x381f81a3820bb853),
                    c1: Fp(0x1cb62c44e2dc8776),
                },
                c2: Fp2 {
                    c0: Fp(0x30657c08bab6d7cf),
                    c1: Fp(0x380f93f2c3012faa),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x30a2e880fed2d926),
                    c1: Fp(0x278634dfaa1f61ca),
                },
                c1: Fp2 {
                    c0: Fp(0x135d0e99d77e5964),
                    c1: Fp(0x2d9d0302674b1ebe),
                },
                c2: Fp2 {
                    c0: Fp(0x3431649aa5919fbf),
                    c1: Fp(0x26e4fb11969ba74f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25ebe1067197d8eb),
                    c1: Fp(0x22a3420446c1b8fe),
                },
                c1: Fp2 {
                    c0: Fp(0x330df8b2986df003),
                    c1: Fp(0x3bde4a3ba9a76d01),
                },
                c2: Fp2 {
                    c0: Fp(0x14cb1ad312d6fafd),
                    c1: Fp(0x23101b670a373b83),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x17c5be9e63c8073b),
                    c1: Fp(0x3d55db029d7b62ec),
                },
                c1: Fp2 {
                    c0: Fp(0x2cb73a99b37c577c),
                    c1: Fp(0x15748f005d292447),
                },
                c2: Fp2 {
                    c0: Fp(0xd3feb4c00f4c582),
                    c1: Fp(0x169f01505f408e34),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3139fb55b887425),
                    c1: Fp(0x222abe45d11d4a83),
                },
                c1: Fp2 {
                    c0: Fp(0x1b11706d5c8402e3),
                    c1: Fp(0x1d74171322ebc4ae),
                },
                c2: Fp2 {
                    c0: Fp(0x62af604f15a7562),
                    c1: Fp(0x1bbd08887fa7fc0a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25c6ed0278ed2155),
                    c1: Fp(0xbccd6638d2d8b7),
                },
                c1: Fp2 {
                    c0: Fp(0x1fa871104b28c341),
                    c1: Fp(0x3d2e43af4463be9b),
                },
                c2: Fp2 {
                    c0: Fp(0xe4c666d16719141),
                    c1: Fp(0x28d2025e5850822a),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3acb03b377c9f7b2),
                    c1: Fp(0x39bafd749c27236f),
                },
                c1: Fp2 {
                    c0: Fp(0xef86f199f16551d),
                    c1: Fp(0x1da8d84a3edfb41d),
                },
                c2: Fp2 {
                    c0: Fp(0x2b5cebbb4aaa5c9a),
                    c1: Fp(0x285f9c6c1a87da90),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x12802c083f2a1bf8),
                    c1: Fp(0x4018cd2e8db7afc9),
                },
                c1: Fp2 {
                    c0: Fp(0x2e84f0f0158c96bc),
                    c1: Fp(0x2bcc174d9bf0860),
                },
                c2: Fp2 {
                    c0: Fp(0x2fd47e7df7eb07c8),
                    c1: Fp(0x3202b433c64e6039),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xf3afe5791296fa2),
                    c1: Fp(0x2fdb2a033711fb72),
                },
                c1: Fp2 {
                    c0: Fp(0x2f411feea96c84d9),
                    c1: Fp(0x166a554d8a606282),
                },
                c2: Fp2 {
                    c0: Fp(0x33ff39f4d0fb9901),
                    c1: Fp(0x396c484909713e00),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2969571449d90c9b),
                    c1: Fp(0x3061b0aafe64a7e9),
                },
                c1: Fp2 {
                    c0: Fp(0xcb7c94fa8f03aa2),
                    c1: Fp(0x65e22943e13d7e0),
                },
                c2: Fp2 {
                    c0: Fp(0x19eb521f3fe1b3a8),
                    c1: Fp(0x2964a729ddd25dc1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x19443ef25f20de39),
                    c1: Fp(0x24a8675fc7352594),
                },
                c1: Fp2 {
                    c0: Fp(0x27b333a0d808348c),
                    c1: Fp(0x11213ea7b5b4199c),
                },
                c2: Fp2 {
                    c0: Fp(0x3cba514d80069a6b),
                    c1: Fp(0x2fa2b04721ad1648),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x11fd94c793f2ce61),
                    c1: Fp(0x3cdf2f2ec2b23b53),
                },
                c1: Fp2 {
                    c0: Fp(0x1cb219281c60729),
                    c1: Fp(0x3561eda7736dd573),
                },
                c2: Fp2 {
                    c0: Fp(0x22c59b5f46ae28ca),
                    c1: Fp(0x25fb032500d20d22),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x278b7b6f4e2cacae),
                    c1: Fp(0x3081f1c598bca01b),
                },
                c1: Fp2 {
                    c0: Fp(0x3d9d0d1d4b3b8267),
                    c1: Fp(0xdc908b4ff96b12d),
                },
                c2: Fp2 {
                    c0: Fp(0x27c86aa886504ce8),
                    c1: Fp(0x23b7c3ec1d6980a0),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ea788ab81a0fd9a),
                    c1: Fp(0x2486c6e534c68c64),
                },
                c1: Fp2 {
                    c0: Fp(0x25cfe310deeb8588),
                    c1: Fp(0xea4c851e3f51702),
                },
                c2: Fp2 {
                    c0: Fp(0x3b4dedae16c62ee9),
                    c1: Fp(0x3afa6e666ea75348),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3868a9327ca33149),
                    c1: Fp(0x103379569ef93d03),
                },
                c1: Fp2 {
                    c0: Fp(0x1726bdfdad6152a),
                    c1: Fp(0x319f6c8e34aeaba2),
                },
                c2: Fp2 {
                    c0: Fp(0x1a2cef24df807760),
                    c1: Fp(0x20b12f2fb2d0e9b3),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3d245cf04f5791b),
                    c1: Fp(0xca64828d6843170),
                },
                c1: Fp2 {
                    c0: Fp(0xaa2a223f9411476),
                    c1: Fp(0x1629f961c1958c99),
                },
                c2: Fp2 {
                    c0: Fp(0x4111b26550884190),
                    c1: Fp(0x38611b4faa31892b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1593fa390af3e41e),
                    c1: Fp(0x1576d62862411632),
                },
                c1: Fp2 {
                    c0: Fp(0x151cd7c58c1efa41),
                    c1: Fp(0x2cf0479c3e878103),
                },
                c2: Fp2 {
                    c0: Fp(0x3835f28bbd67bfb7),
                    c1: Fp(0x151e7219bc402465),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x106afb8cbfcb1c19),
                    c1: Fp(0x4024c483b55f7a9d),
                },
                c1: Fp2 {
                    c0: Fp(0x34beffdb82238aca),
                    c1: Fp(0x3349f3d3346d3eb),
                },
                c2: Fp2 {
                    c0: Fp(0x3c90bad089e45b5a),
                    c1: Fp(0x30c6e8e0676397f0),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2ecdd4244cab2682),
                    c1: Fp(0x26ccebe43d02fa08),
                },
                c1: Fp2 {
                    c0: Fp(0x1bdf568b5abd822e),
                    c1: Fp(0x1c72fec1983c24b9),
                },
                c2: Fp2 {
                    c0: Fp(0x32cde7da5a2e7c04),
                    c1: Fp(0x35ed7b9b5444f32c),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x244ecf2565296a73),
                    c1: Fp(0x158ef4064604bd0e),
                },
                c1: Fp2 {
                    c0: Fp(0x161a7d5be454922e),
                    c1: Fp(0x351d7d3f45ec7bd6),
                },
                c2: Fp2 {
                    c0: Fp(0x2e8fe708a1bcd8c3),
                    c1: Fp(0x213c3ba3e9f1dfa2),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x371ef979730f229c),
                    c1: Fp(0x3244140f23f57d2a),
                },
                c1: Fp2 {
                    c0: Fp(0x27725837096104c0),
                    c1: Fp(0x10cafa232a61d8c5),
                },
                c2: Fp2 {
                    c0: Fp(0x1ce40b6fbbf88b54),
                    c1: Fp(0xc71abed0f3429e4),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x187b027472b89360),
                    c1: Fp(0x1e87db21f34a946f),
                },
                c1: Fp2 {
                    c0: Fp(0x1a2dfb40ea68f923),
                    c1: Fp(0x24e21b8dc660feea),
                },
                c2: Fp2 {
                    c0: Fp(0xe2b721e3f0f6a0b),
                    c1: Fp(0x1501bee4f2739f64),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2c190831060d2c3f),
                    c1: Fp(0x8106685876814ca),
                },
                c1: Fp2 {
                    c0: Fp(0x13008c33bdf43e64),
                    c1: Fp(0x34e397fe69bdef76),
                },
                c2: Fp2 {
                    c0: Fp(0x17c25f35a10bacc0),
                    c1: Fp(0x2d2f85e4610580fa),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x344e347c90910a35),
                    c1: Fp(0x270b3bb9ae08cdd9),
                },
                c1: Fp2 {
                    c0: Fp(0x2e9afa51eca61db4),
                    c1: Fp(0x33f3fb71058e55ec),
                },
                c2: Fp2 {
                    c0: Fp(0x134dd53bdeffd836),
                    c1: Fp(0x3c22c1c5f94a82e4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x82c95604367dd63),
                    c1: Fp(0x18eb43838bffd616),
                },
                c1: Fp2 {
                    c0: Fp(0xa43fa9b277a9201),
                    c1: Fp(0x343bfbc46372bd06),
                },
                c2: Fp2 {
                    c0: Fp(0x1873b7c40f058862),
                    c1: Fp(0x248cc6d14fe66c4d),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x180d3c6028d3c445),
                    c1: Fp(0x337d79cadde79058),
                },
                c1: Fp2 {
                    c0: Fp(0x2e4e08874cfbf8cd),
                    c1: Fp(0x190f8ec01b4e3994),
                },
                c2: Fp2 {
                    c0: Fp(0x33a06241b402850d),
                    c1: Fp(0xe0552d2809f2478),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x302e96910e19a8e7),
                    c1: Fp(0xb9ee0f3b65cd639),
                },
                c1: Fp2 {
                    c0: Fp(0xa4e8355e5e0fde5),
                    c1: Fp(0x2f21d6c7669e9d10),
                },
                c2: Fp2 {
                    c0: Fp(0x119ebe7bad6dcdb7),
                    c1: Fp(0x2dcc54ab42cbcf19),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x214142d066497e3b),
                    c1: Fp(0x26fc29b5e269fbd5),
                },
                c1: Fp2 {
                    c0: Fp(0xcab7d6fdd79b53f),
                    c1: Fp(0x47665ad8c0dbe41),
                },
                c2: Fp2 {
                    c0: Fp(0x1dc0dc89de799be0),
                    c1: Fp(0x2a380f311e216d4f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x313ecc609b0bba65),
                    c1: Fp(0x26034fdde1e4bb00),
                },
                c1: Fp2 {
                    c0: Fp(0x277e614cd68f1577),
                    c1: Fp(0x12165771786fda0d),
                },
                c2: Fp2 {
                    c0: Fp(0x2b21d78a77a25e89),
                    c1: Fp(0xe911edb7cc9aa4c),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xa66868701499899),
                    c1: Fp(0xcdbe76635826edb),
                },
                c1: Fp2 {
                    c0: Fp(0xbd96d0df4a9869c),
                    c1: Fp(0xf3249ce7c748f45),
                },
                c2: Fp2 {
                    c0: Fp(0x31447f390ca58f40),
                    c1: Fp(0x168c7bfd69f6aae8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x29ea28b000c1216a),
                    c1: Fp(0xc006a512ed86383),
                },
                c1: Fp2 {
                    c0: Fp(0x2184bfd8eba79752),
                    c1: Fp(0x3c30ec0a692a22a3),
                },
                c2: Fp2 {
                    c0: Fp(0x2a87baec268ea488),
                    c1: Fp(0x36305263ab7cc76e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x22914d5bc6aa3426),
                    c1: Fp(0x3c14c4f5e0e1890c),
                },
                c1: Fp2 {
                    c0: Fp(0x36e7ae8e3fdcbfda),
                    c1: Fp(0x3947687ae16cab),
                },
                c2: Fp2 {
                    c0: Fp(0x376400f381955e39),
                    c1: Fp(0x3053a1537900f436),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3b59e828e9287248),
                    c1: Fp(0x31e3b418c03f0e08),
                },
                c1: Fp2 {
                    c0: Fp(0x45dbead6ea5fc18),
                    c1: Fp(0x326413511203a603),
                },
                c2: Fp2 {
                    c0: Fp(0x23786732e5009717),
                    c1: Fp(0x26c4630d631cd64),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25bcf58b363cc51d),
                    c1: Fp(0x17a4f0b12267dc85),
                },
                c1: Fp2 {
                    c0: Fp(0x20b2cb3a0dc737b5),
                    c1: Fp(0x3d7c34c547acae60),
                },
                c2: Fp2 {
                    c0: Fp(0x3e524f5a78ee7474),
                    c1: Fp(0x27002435a4294378),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3eb8ae98e91b8f20),
                    c1: Fp(0x733e167672758f4),
                },
                c1: Fp2 {
                    c0: Fp(0x4014e8ef326e4a58),
                    c1: Fp(0x333085981d30d8fc),
                },
                c2: Fp2 {
                    c0: Fp(0x2f491e6e7154cdb),
                    c1: Fp(0xf85425f5d9eb038),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1f4f97c2d77c8bfb),
                    c1: Fp(0x2e7e902971094bd7),
                },
                c1: Fp2 {
                    c0: Fp(0xff8a6dfb6a55a9),
                    c1: Fp(0x272746de279138cc),
                },
                c2: Fp2 {
                    c0: Fp(0x318e2d12f350a344),
                    c1: Fp(0x273fe5461026d21c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x34d3f9ab54d1672f),
                    c1: Fp(0x529b528fea7b0f),
                },
                c1: Fp2 {
                    c0: Fp(0x3992bd2d65a12253),
                    c1: Fp(0x90e3d46bebee8db),
                },
                c2: Fp2 {
                    c0: Fp(0x3c3f3714ca7a97df),
                    c1: Fp(0xbd4101bf4443527),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x401631cb86625526),
                    c1: Fp(0x3f35e77f7346d156),
                },
                c1: Fp2 {
                    c0: Fp(0x40dd1ef63b601150),
                    c1: Fp(0xae6344e5ec1a073),
                },
                c2: Fp2 {
                    c0: Fp(0x18a70e4880ffc5cb),
                    c1: Fp(0x1dd84414961aa815),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xac1acdfb54236fc),
                    c1: Fp(0xd0c2c4641313285),
                },
                c1: Fp2 {
                    c0: Fp(0x3d44ee3951e1ec5e),
                    c1: Fp(0xb1836b4b5e051a5),
                },
                c2: Fp2 {
                    c0: Fp(0x689844c1f190793),
                    c1: Fp(0x6ee39d695056ba8),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c5f91aeac99a0a9),
                    c1: Fp(0x3214a54b7831e81b),
                },
                c1: Fp2 {
                    c0: Fp(0x2ec352d496490f22),
                    c1: Fp(0xadfaff9b49c0d96),
                },
                c2: Fp2 {
                    c0: Fp(0x3f1fa3fc5b207db0),
                    c1: Fp(0x215783e03ce9aab3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2a047cb9da6c2bae),
                    c1: Fp(0x3f2b095404425af5),
                },
                c1: Fp2 {
                    c0: Fp(0x581cb6502c4db00),
                    c1: Fp(0x3c05eaf0085be05),
                },
                c2: Fp2 {
                    c0: Fp(0x13ced73b74494c2e),
                    c1: Fp(0x180261e3a463302a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x169673bad78ae947),
                    c1: Fp(0x34e29fa8c2cf81b3),
                },
                c1: Fp2 {
                    c0: Fp(0x2cfa21135708e258),
                    c1: Fp(0x3414c36fabf16e52),
                },
                c2: Fp2 {
                    c0: Fp(0x119b953c2c747b57),
                    c1: Fp(0x3f9f9f38c4d29bf0),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x55eeb99683cf25a),
                    c1: Fp(0x3cff0192378e5edc),
                },
                c1: Fp2 {
                    c0: Fp(0x2938548504f59da4),
                    c1: Fp(0x4099ffbdcc636602),
                },
                c2: Fp2 {
                    c0: Fp(0x373858cb0a04df2b),
                    c1: Fp(0x29ce781a113b6743),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1c8ddffd95b1bc32),
                    c1: Fp(0x247790d1c72afd7f),
                },
                c1: Fp2 {
                    c0: Fp(0x168609c273cd07e),
                    c1: Fp(0x1a3e56ab9cc63544),
                },
                c2: Fp2 {
                    c0: Fp(0x14430d346f512466),
                    c1: Fp(0x30a719a27b5b7c9f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a9ec2411bfc7c3),
                    c1: Fp(0x30860ae318332f36),
                },
                c1: Fp2 {
                    c0: Fp(0x15f01698b82efeef),
                    c1: Fp(0x360e09a61e03a448),
                },
                c2: Fp2 {
                    c0: Fp(0x1b9667ef3680ed37),
                    c1: Fp(0x15feafe937400328),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x215936fc162af248),
                    c1: Fp(0xfde1d603524283f),
                },
                c1: Fp2 {
                    c0: Fp(0xd1548322d55c54b),
                    c1: Fp(0x166d1fd15fd14cde),
                },
                c2: Fp2 {
                    c0: Fp(0x3df68b7b9275eabc),
                    c1: Fp(0x2fe582bfffdba90c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x18f349654095ed73),
                    c1: Fp(0x4032d1739688ecb7),
                },
                c1: Fp2 {
                    c0: Fp(0x39007dc086ec5f4e),
                    c1: Fp(0x1972c2329673a472),
                },
                c2: Fp2 {
                    c0: Fp(0x372dc5d44acf6ae3),
                    c1: Fp(0x3c6f7aad7b94d5af),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1b13d4f2f98ac7f5),
                    c1: Fp(0x1880bdfbdf5d0c10),
                },
                c1: Fp2 {
                    c0: Fp(0x24944c895ed728ba),
                    c1: Fp(0x2fd800917d6cda02),
                },
                c2: Fp2 {
                    c0: Fp(0x116953de58e2c248),
                    c1: Fp(0xf69d43674298d47),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2164a93b9b2b9756),
                    c1: Fp(0x6609f55ac4fda6d),
                },
                c1: Fp2 {
                    c0: Fp(0x1e08df7203f79888),
                    c1: Fp(0x3f26b11273b99062),
                },
                c2: Fp2 {
                    c0: Fp(0x195dcc6556d24e0d),
                    c1: Fp(0x3cd4da7d15eb1481),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xa089d756ac1e984),
                    c1: Fp(0x1296166bc7de70f1),
                },
                c1: Fp2 {
                    c0: Fp(0x35b8f894985162a7),
                    c1: Fp(0x3fcc60a5bc2b9abb),
                },
                c2: Fp2 {
                    c0: Fp(0x184b88c71273506b),
                    c1: Fp(0xf466bd19731cb01),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x20d7fcadcf10afc2),
                    c1: Fp(0x17a929a436d6fd97),
                },
                c1: Fp2 {
                    c0: Fp(0x2606ebbfec3f1dd5),
                    c1: Fp(0x2dea67e3413d8ccb),
                },
                c2: Fp2 {
                    c0: Fp(0x2facca6ca33859cf),
                    c1: Fp(0x2aa6db0b236fdf9a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x16c4dc34036f3d92),
                    c1: Fp(0x50028210422445a),
                },
                c1: Fp2 {
                    c0: Fp(0x3658698e08131f3a),
                    c1: Fp(0x32b7b0b22da26e4a),
                },
                c2: Fp2 {
                    c0: Fp(0x7518a51ef5150d6),
                    c1: Fp(0x2e148576a70ce532),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x83e6ea38842bde),
                    c1: Fp(0x37fb700e2e543e30),
                },
                c1: Fp2 {
                    c0: Fp(0x3573f69b98efd97f),
                    c1: Fp(0x400e5eab4cf9df99),
                },
                c2: Fp2 {
                    c0: Fp(0x3e6c48da9bf3fe22),
                    c1: Fp(0x224b477e03e1ff67),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3764b3d5057cf31a),
                    c1: Fp(0x1bfb07f58ce30d88),
                },
                c1: Fp2 {
                    c0: Fp(0x3d655450c9072140),
                    c1: Fp(0x31d01686a8116541),
                },
                c2: Fp2 {
                    c0: Fp(0x1e6fa37916d6a937),
                    c1: Fp(0x1a0b684fd97dc25d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x27355dab59df4b08),
                    c1: Fp(0x2c87ee0e0d2dc679),
                },
                c1: Fp2 {
                    c0: Fp(0x270cc0b747b50457),
                    c1: Fp(0x344f65ca206f6882),
                },
                c2: Fp2 {
                    c0: Fp(0x2d13311f4bf170ca),
                    c1: Fp(0x24dc30c0b10f18d5),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x272ec1b9a902ee2d),
                    c1: Fp(0x2278e9188ea77549),
                },
                c1: Fp2 {
                    c0: Fp(0xd9c843a4c95eb98),
                    c1: Fp(0x5219c5bdb791b64),
                },
                c2: Fp2 {
                    c0: Fp(0x1eeecadfaacd9862),
                    c1: Fp(0x1203c514ed57c598),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x11b443673d52740),
                    c1: Fp(0x1eb2db54542dbae9),
                },
                c1: Fp2 {
                    c0: Fp(0x18b69e69fa940067),
                    c1: Fp(0x2bfd127acd768fc0),
                },
                c2: Fp2 {
                    c0: Fp(0x1d6b031f72209960),
                    c1: Fp(0x3ddd6cd1b4e79338),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x389f7095e1aad15b),
                    c1: Fp(0x2c6d03d72461fb78),
                },
                c1: Fp2 {
                    c0: Fp(0x2f6209f0561b5984),
                    c1: Fp(0x2668ee598e92cddf),
                },
                c2: Fp2 {
                    c0: Fp(0x3326b3cae179b4b1),
                    c1: Fp(0x40eddd2a589f1cb),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xf4bedb7c3c2decc),
                    c1: Fp(0x21913aea9db4710f),
                },
                c1: Fp2 {
                    c0: Fp(0xac95a18e7d19050),
                    c1: Fp(0x12279a6dcb559388),
                },
                c2: Fp2 {
                    c0: Fp(0x6a02e50aec95825),
                    c1: Fp(0x2b6f750f181cc9c4),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1f9a5a2984b28a8b),
                    c1: Fp(0x84951699ca9238a),
                },
                c1: Fp2 {
                    c0: Fp(0x1621c77c3676993e),
                    c1: Fp(0x8265af812f52bb5),
                },
                c2: Fp2 {
                    c0: Fp(0x3e77ff6d436c2000),
                    c1: Fp(0xc25f2d715a8a7d1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1d2a37c9c7e872fe),
                    c1: Fp(0x6856278a6aac9b1),
                },
                c1: Fp2 {
                    c0: Fp(0x39bc4eaf2d3ef2d0),
                    c1: Fp(0x24212615bf483525),
                },
                c2: Fp2 {
                    c0: Fp(0x1e9550764e3bea7f),
                    c1: Fp(0x2960983768ab241b),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1a43732400613507),
                    c1: Fp(0xcc0f34a1fe93b82),
                },
                c1: Fp2 {
                    c0: Fp(0x8fe7d2be82bfe3e),
                    c1: Fp(0x3aa33198fd633294),
                },
                c2: Fp2 {
                    c0: Fp(0x3cace4aa236e78ab),
                    c1: Fp(0x3d5c8bd326d587ea),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ad13b6d4d34fec7),
                    c1: Fp(0x122c9ba548430cef),
                },
                c1: Fp2 {
                    c0: Fp(0xe801b5457a9c6db),
                    c1: Fp(0x1bdfdd822c654067),
                },
                c2: Fp2 {
                    c0: Fp(0x3143ed347ee3a0ef),
                    c1: Fp(0x40404301317e6761),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x6c580b97a8450c2),
                    c1: Fp(0x6c2ba97e8334571),
                },
                c1: Fp2 {
                    c0: Fp(0x352ed748699e4ca8),
                    c1: Fp(0xfa20743fe3220f5),
                },
                c2: Fp2 {
                    c0: Fp(0x25e5d48127929cbe),
                    c1: Fp(0x413482b2555455c4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x15f095df0f732b60),
                    c1: Fp(0x3a645787102d3b35),
                },
                c1: Fp2 {
                    c0: Fp(0x23a0ba1248d958bb),
                    c1: Fp(0x19a906dae526d66d),
                },
                c2: Fp2 {
                    c0: Fp(0x1543f5ee2dad3ac5),
                    c1: Fp(0x3dc5e80b9d606f72),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xaf9a0df859d515e),
                    c1: Fp(0x3889dd1a9df80df2),
                },
                c1: Fp2 {
                    c0: Fp(0xbe60105590f619c),
                    c1: Fp(0x660a30327cced49),
                },
                c2: Fp2 {
                    c0: Fp(0x33303c9eb74db1fe),
                    c1: Fp(0x271bf62c44e09078),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x397bf3427f6a4d1d),
                    c1: Fp(0x1b61ba447be16d1c),
                },
                c1: Fp2 {
                    c0: Fp(0x72cb59a7b8c85ee),
                    c1: Fp(0xae57135c8ab29a5),
                },
                c2: Fp2 {
                    c0: Fp(0x40b8e323229c01e5),
                    c1: Fp(0x29ad49fd22d3af0f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a30a1ed511eb68e),
                    c1: Fp(0xdcfe0d53243ac30),
                },
                c1: Fp2 {
                    c0: Fp(0xac8e3c3a44cf165),
                    c1: Fp(0x2d3b90f6bf5b784e),
                },
                c2: Fp2 {
                    c0: Fp(0x1029db0d1f27df5a),
                    c1: Fp(0x2ade47b38af9b788),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x18a589d4a1be3abe),
                    c1: Fp(0x32017c3618ddb2ef),
                },
                c1: Fp2 {
                    c0: Fp(0x7c86dc61f252c62),
                    c1: Fp(0x11844b931884a569),
                },
                c2: Fp2 {
                    c0: Fp(0x3f8b04896bffd3da),
                    c1: Fp(0xc1d9992351b7592),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1694dfd3994645e7),
                    c1: Fp(0x93d33acd196be04),
                },
                c1: Fp2 {
                    c0: Fp(0x24173b73bf1c3cc4),
                    c1: Fp(0xc6f164adb2b2888),
                },
                c2: Fp2 {
                    c0: Fp(0x152f72a98d334688),
                    c1: Fp(0x34a48b480bf4df3a),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2ae5466abb6afb09),
                    c1: Fp(0x1549261654b700c4),
                },
                c1: Fp2 {
                    c0: Fp(0x185a59055cf6c363),
                    c1: Fp(0x274a5f496b79e32),
                },
                c2: Fp2 {
                    c0: Fp(0x1bf454a0898c0a53),
                    c1: Fp(0x3820162a17adf7ec),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x4bba632854b23cb),
                    c1: Fp(0x120d53b7ad0c75e3),
                },
                c1: Fp2 {
                    c0: Fp(0x2ddac1001a993ed5),
                    c1: Fp(0xd9df73e9167f18a),
                },
                c2: Fp2 {
                    c0: Fp(0x11b6e366c1429d8b),
                    c1: Fp(0xcd9e78d02cd8a74),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1a580972bfcea685),
                    c1: Fp(0x3b7632ca94a18585),
                },
                c1: Fp2 {
                    c0: Fp(0x1d1f6e6ca97de9d4),
                    c1: Fp(0x28c3e56670df3f84),
                },
                c2: Fp2 {
                    c0: Fp(0x325de370c7d7bf77),
                    c1: Fp(0x375ad5f75eba4027),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x36b3b7cb5fc90e46),
                    c1: Fp(0x2e94e81a39295e6c),
                },
                c1: Fp2 {
                    c0: Fp(0x35953e77b914c885),
                    c1: Fp(0x298fdf37abe95279),
                },
                c2: Fp2 {
                    c0: Fp(0x19c619ee1fef72e9),
                    c1: Fp(0x16413dfb36073ee),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x4df2ef1e3f7d288),
                    c1: Fp(0x358747af9a984e0c),
                },
                c1: Fp2 {
                    c0: Fp(0x1ef48af0a74ce64d),
                    c1: Fp(0xda8bc31d9232d6e),
                },
                c2: Fp2 {
                    c0: Fp(0x1ed307571753bbc2),
                    c1: Fp(0x2f98aa6812d7543d),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1644afb018c599aa),
                    c1: Fp(0x3207cf7c19c274ce),
                },
                c1: Fp2 {
                    c0: Fp(0x2314db4152872d19),
                    c1: Fp(0x21b68cbf3482eb15),
                },
                c2: Fp2 {
                    c0: Fp(0x3459e6421f2b6022),
                    c1: Fp(0x5421a2fce2a125e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd938bea5ec8c008),
                    c1: Fp(0x2ae7f08515164116),
                },
                c1: Fp2 {
                    c0: Fp(0x1f6d57842ffc829d),
                    c1: Fp(0x232b6e20ac1c19b8),
                },
                c2: Fp2 {
                    c0: Fp(0x1c484eca596379e1),
                    c1: Fp(0x5b2060e34b799a3),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xc5270f2d698c3d4),
                    c1: Fp(0x3fd12d26cc8a78ed),
                },
                c1: Fp2 {
                    c0: Fp(0x230241d7c0a3d990),
                    c1: Fp(0x1865630304ca1b09),
                },
                c2: Fp2 {
                    c0: Fp(0x23a9a7a68e34f011),
                    c1: Fp(0x1c691c361005b6a3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe908a245ebbb4a2),
                    c1: Fp(0xc08f5543b35125),
                },
                c1: Fp2 {
                    c0: Fp(0x1609f1b7e5d29e42),
                    c1: Fp(0x2ec31802d00c22b0),
                },
                c2: Fp2 {
                    c0: Fp(0x1e3f563e8f237fe1),
                    c1: Fp(0x2f69ce8f8a6e2930),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x175d228cc8791861),
                    c1: Fp(0x3125fdfb7488c7a5),
                },
                c1: Fp2 {
                    c0: Fp(0x36aba630a6e5e882),
                    c1: Fp(0x40db5b8da5b85cca),
                },
                c2: Fp2 {
                    c0: Fp(0x3c12ed5f86e01395),
                    c1: Fp(0x153672f528e6bfa8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x33804615dc7caac5),
                    c1: Fp(0x164156b0ddc49d0c),
                },
                c1: Fp2 {
                    c0: Fp(0xbc27d458b9ba11),
                    c1: Fp(0x31117d7ba1c26055),
                },
                c2: Fp2 {
                    c0: Fp(0x338583504a968eba),
                    c1: Fp(0x5147f52e2b95ea7),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2437039c73e9e34),
                    c1: Fp(0x17486de6f4a6f04c),
                },
                c1: Fp2 {
                    c0: Fp(0x791fb1d4d40e127),
                    c1: Fp(0x13d5cc085e16f534),
                },
                c2: Fp2 {
                    c0: Fp(0x2c362376c24a268a),
                    c1: Fp(0x3cd99a3a83be6fd7),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2f846101f7b8f64a),
                    c1: Fp(0x2eee35f84cf0dee9),
                },
                c1: Fp2 {
                    c0: Fp(0x3e21b202d697e7b5),
                    c1: Fp(0x267fbc158e1cec26),
                },
                c2: Fp2 {
                    c0: Fp(0xe3266d36eb8ac13),
                    c1: Fp(0x3103b57c8e678863),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xbe96d0a9ca488),
                    c1: Fp(0x26d61a5ae7a5a72e),
                },
                c1: Fp2 {
                    c0: Fp(0x26d8868a5a7ef341),
                    c1: Fp(0xdd50fbf9c83efbe),
                },
                c2: Fp2 {
                    c0: Fp(0x1749476f32031f96),
                    c1: Fp(0x3062c477e5f66b9),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x35ad30da5df64b72),
                    c1: Fp(0x31c2fee900fb2532),
                },
                c1: Fp2 {
                    c0: Fp(0x222ec93742162907),
                    c1: Fp(0x11edf5385b612a2b),
                },
                c2: Fp2 {
                    c0: Fp(0x5b1b72c71f87cda),
                    c1: Fp(0xd0d4d6130a7a3fc),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2565bccd962a373e),
                    c1: Fp(0x34fb394e6b5f3a71),
                },
                c1: Fp2 {
                    c0: Fp(0x5a0b3121824731b),
                    c1: Fp(0x12db942ad516d8fa),
                },
                c2: Fp2 {
                    c0: Fp(0x22d4361213cab5f5),
                    c1: Fp(0x3f485b98d79b8cf8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd7071fcccca1c55),
                    c1: Fp(0x1ea16094725652ff),
                },
                c1: Fp2 {
                    c0: Fp(0xcaeaa993b624e01),
                    c1: Fp(0x2ab0faf00bd000df),
                },
                c2: Fp2 {
                    c0: Fp(0x402f7c104abc7dee),
                    c1: Fp(0x26ea588451a649b7),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25ef206a4fd3b7fc),
                    c1: Fp(0x8637b9753c1d2c1),
                },
                c1: Fp2 {
                    c0: Fp(0x10fd8054cb4633ae),
                    c1: Fp(0x188af512d13bd44),
                },
                c2: Fp2 {
                    c0: Fp(0x13e7d143b7087b17),
                    c1: Fp(0x2e061b8a5c8942ac),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x393f71d1694e1781),
                    c1: Fp(0x2157f45c96e7ec28),
                },
                c1: Fp2 {
                    c0: Fp(0x3d3be5abfe5e1b51),
                    c1: Fp(0x22d5168f1ff5a381),
                },
                c2: Fp2 {
                    c0: Fp(0x2b0249ae2a7fd11c),
                    c1: Fp(0x54d4c272a29e460),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2331893951cb9642),
                    c1: Fp(0x37b4d5f69e72257e),
                },
                c1: Fp2 {
                    c0: Fp(0x3b5fd835dcb69d00),
                    c1: Fp(0x415613f33c0ef58b),
                },
                c2: Fp2 {
                    c0: Fp(0x277cca4f2fb38003),
                    c1: Fp(0x4041d965eadcbf3b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c9231304e501530),
                    c1: Fp(0x3295d712df85b9c7),
                },
                c1: Fp2 {
                    c0: Fp(0x1a795485a175c67e),
                    c1: Fp(0x23dab52a2b2006f6),
                },
                c2: Fp2 {
                    c0: Fp(0x1f1f3fd18bbbd0e9),
                    c1: Fp(0x2f3bc11e4b8a9380),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x29296bdcbb66a0b6),
                    c1: Fp(0x33a94ea9ae783b6e),
                },
                c1: Fp2 {
                    c0: Fp(0x35ef178dba68adcd),
                    c1: Fp(0x3f22aa3cc38dc6b5),
                },
                c2: Fp2 {
                    c0: Fp(0x175f17f3e9e071ed),
                    c1: Fp(0x31be24737d979804),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x266819086d7be65e),
                    c1: Fp(0x16685f5401b386cf),
                },
                c1: Fp2 {
                    c0: Fp(0x2a31b2d8a18e2abf),
                    c1: Fp(0x10dbadfa223d3eba),
                },
                c2: Fp2 {
                    c0: Fp(0x326885393f7fe264),
                    c1: Fp(0x2de3e729ea069f2a),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x109f504214cb9720),
                    c1: Fp(0x1ac84f5cda01f39b),
                },
                c1: Fp2 {
                    c0: Fp(0x14ef81ad95f81cb3),
                    c1: Fp(0x6d07e5e1757794a),
                },
                c2: Fp2 {
                    c0: Fp(0x3b3cc1a62132a4ee),
                    c1: Fp(0x7bdcd76287820c8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x37488af85a9ea00e),
                    c1: Fp(0x17c2a2efeb621d6e),
                },
                c1: Fp2 {
                    c0: Fp(0x2e241487e65e27ea),
                    c1: Fp(0x14d5527520b656eb),
                },
                c2: Fp2 {
                    c0: Fp(0x368bda2a5c1a66c3),
                    c1: Fp(0x1ec4275f3306815a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2d9690b0658388cd),
                    c1: Fp(0x33fc4d0b5154ac81),
                },
                c1: Fp2 {
                    c0: Fp(0x1c1d882f07715ec9),
                    c1: Fp(0x29bac75a893ce5f0),
                },
                c2: Fp2 {
                    c0: Fp(0x341ee936133b0e11),
                    c1: Fp(0x4d4661575fe9cb4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x8ef65afb3496bd6),
                    c1: Fp(0x289ebc6d9cff5c5e),
                },
                c1: Fp2 {
                    c0: Fp(0x3d41f02cae7f29ee),
                    c1: Fp(0x154d06051a8d0b3c),
                },
                c2: Fp2 {
                    c0: Fp(0xe75555e33f0c3ef),
                    c1: Fp(0x613bc5214a765ef),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x18e4cb839962a979),
                    c1: Fp(0x3f575b5a11f83c01),
                },
                c1: Fp2 {
                    c0: Fp(0x3bfea95a6be5d150),
                    c1: Fp(0x18d702458476d07f),
                },
                c2: Fp2 {
                    c0: Fp(0x414d63537e9f903f),
                    c1: Fp(0x7f55773cfabea73),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x22034ab3b846641a),
                    c1: Fp(0x3b0f05649af2f0d4),
                },
                c1: Fp2 {
                    c0: Fp(0x3ab691df156028d3),
                    c1: Fp(0x26e19266976cdba9),
                },
                c2: Fp2 {
                    c0: Fp(0x3adaf32a324f4ead),
                    c1: Fp(0x2397fadefcd8cb0c),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3d04a2108f48a6cc),
                    c1: Fp(0x1b49d2199144f55c),
                },
                c1: Fp2 {
                    c0: Fp(0x3d82d59627407483),
                    c1: Fp(0x9631f621d676ecf),
                },
                c2: Fp2 {
                    c0: Fp(0xf6fe7c87873c084),
                    c1: Fp(0x3d14bf427604180e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1bee7b17955e6390),
                    c1: Fp(0x40dae91ffb85c983),
                },
                c1: Fp2 {
                    c0: Fp(0x2af3c0b675d1293),
                    c1: Fp(0x3347d02dd2b4ce6d),
                },
                c2: Fp2 {
                    c0: Fp(0x404cb3773637141b),
                    c1: Fp(0xe846f7e7da3ed97),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3db669157e3908f4),
                    c1: Fp(0x1eea1d29d4539ed6),
                },
                c1: Fp2 {
                    c0: Fp(0xd3a3e88d8635237),
                    c1: Fp(0x6ecc85502499f77),
                },
                c2: Fp2 {
                    c0: Fp(0x2ca1055e587ecc3a),
                    c1: Fp(0x2265ea1ba7bfb975),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xecb971819d38731),
                    c1: Fp(0x30f7825a7d8b2753),
                },
                c1: Fp2 {
                    c0: Fp(0x1cdda9b4c25665a7),
                    c1: Fp(0x406148d6eccf2d93),
                },
                c2: Fp2 {
                    c0: Fp(0x1848b66c8802a14b),
                    c1: Fp(0x28f6d431b3e0bf7e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ae9b8c4d92cfaad),
                    c1: Fp(0x628ed5c0ff9cfb2),
                },
                c1: Fp2 {
                    c0: Fp(0x412c6e2acf512412),
                    c1: Fp(0x1d4218d2d1ef816),
                },
                c2: Fp2 {
                    c0: Fp(0x2ae804c03dd0288),
                    c1: Fp(0xe3e57d6b1c8aa2a),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x10754d8f72465d2a),
                    c1: Fp(0x1e423b55c9780b18),
                },
                c1: Fp2 {
                    c0: Fp(0x1bbfd8ab83f03d49),
                    c1: Fp(0x100613b6b8441431),
                },
                c2: Fp2 {
                    c0: Fp(0x176e7b5587bc1852),
                    c1: Fp(0x20e18d920d1cb16f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x307bd85750329ba1),
                    c1: Fp(0x25ed06add2e29c69),
                },
                c1: Fp2 {
                    c0: Fp(0x9933b0c9d4ab45d),
                    c1: Fp(0xbebfd0781defe1c),
                },
                c2: Fp2 {
                    c0: Fp(0x1fb156db2d386991),
                    c1: Fp(0x16ec1aa4d3c29de3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe9aef29b9125ff7),
                    c1: Fp(0x16a1709b90f19cf8),
                },
                c1: Fp2 {
                    c0: Fp(0x4153f7c138bdd6d8),
                    c1: Fp(0x3c460a29bf69bb43),
                },
                c2: Fp2 {
                    c0: Fp(0x413d899acaa9b43a),
                    c1: Fp(0x1a4b7cc13edc8f1d),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xdd01917f6f4421b),
                    c1: Fp(0x24a6c113628017ea),
                },
                c1: Fp2 {
                    c0: Fp(0x1517cee72d4cd657),
                    c1: Fp(0xa3ed9d81a7a878),
                },
                c2: Fp2 {
                    c0: Fp(0x368bd499c0962972),
                    c1: Fp(0xdfbae5ed1b09a13),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x135a26a7d71cd961),
                    c1: Fp(0x2b27ea521c7709fd),
                },
                c1: Fp2 {
                    c0: Fp(0x361bd0ee502d2323),
                    c1: Fp(0x29537e1e6fbc0901),
                },
                c2: Fp2 {
                    c0: Fp(0x17fd532c85ea9de3),
                    c1: Fp(0x35c69fae21cc5063),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x19ec2111b40a280f),
                    c1: Fp(0x8b16136b3f1c67f),
                },
                c1: Fp2 {
                    c0: Fp(0x1cd1c1472c8c0747),
                    c1: Fp(0x15a000c4560764ca),
                },
                c2: Fp2 {
                    c0: Fp(0x81a8478778e9ffe),
                    c1: Fp(0x400801b350e3851f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ecb094a1b11dbd),
                    c1: Fp(0x39de3c933272790d),
                },
                c1: Fp2 {
                    c0: Fp(0x311fa43910f00d11),
                    c1: Fp(0x19d63201b05fb0d3),
                },
                c2: Fp2 {
                    c0: Fp(0x3b7269b23437529b),
                    c1: Fp(0x1d19e35741458e48),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x202dbe228b7823cc),
                    c1: Fp(0x4369fd0c1c2d62c),
                },
                c1: Fp2 {
                    c0: Fp(0x34796afca2c71a66),
                    c1: Fp(0x7a9de22ce03e42a),
                },
                c2: Fp2 {
                    c0: Fp(0x22e56df42b6103a0),
                    c1: Fp(0xea0bbd9e7816437),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x248e1f152bc4d99),
                    c1: Fp(0xe8b07ad20e2e9cb),
                },
                c1: Fp2 {
                    c0: Fp(0x33f1162bef8abdc3),
                    c1: Fp(0x2410ba2ca82eb0fa),
                },
                c2: Fp2 {
                    c0: Fp(0x1b1205a4c828e2ef),
                    c1: Fp(0x10aed5c0532c471c),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x7918e65c7b64ce0),
                    c1: Fp(0x342c5bf573bf88ed),
                },
                c1: Fp2 {
                    c0: Fp(0x16a4177ec2b14233),
                    c1: Fp(0x326303858345616),
                },
                c2: Fp2 {
                    c0: Fp(0x174f9cd7350cfcdd),
                    c1: Fp(0x274c936439cafea7),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3dc5f4173360f184),
                    c1: Fp(0xec6f669c97a6a2a),
                },
                c1: Fp2 {
                    c0: Fp(0x31f6a044de33a561),
                    c1: Fp(0x324423b92f5be8df),
                },
                c2: Fp2 {
                    c0: Fp(0x35f847845dce55b5),
                    c1: Fp(0x2015987e589f25e4),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xcb21756fa2e2435),
                    c1: Fp(0x1b1463fc3746cb77),
                },
                c1: Fp2 {
                    c0: Fp(0x2e1398770097bec),
                    c1: Fp(0x2cf055dbb4b14282),
                },
                c2: Fp2 {
                    c0: Fp(0x8c8c88acfde56a7),
                    c1: Fp(0xee7465dec5375d1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3975331a4dc4d78),
                    c1: Fp(0x27198d8b12d562a9),
                },
                c1: Fp2 {
                    c0: Fp(0x1bff85b9776a045e),
                    c1: Fp(0x36f15edacaf7bcf1),
                },
                c2: Fp2 {
                    c0: Fp(0x22764663abe6524b),
                    c1: Fp(0x3f94df971ed424f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3bdfdb74ff10fdec),
                    c1: Fp(0x31fcfd99ef623c5e),
                },
                c1: Fp2 {
                    c0: Fp(0x3a624f4afb9e86c1),
                    c1: Fp(0x20037d53c6b41534),
                },
                c2: Fp2 {
                    c0: Fp(0x1f08c2ffc997f876),
                    c1: Fp(0xbd1489bd44cada3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x215cc664fd11ef65),
                    c1: Fp(0x2e0a4a2cb5aa0e79),
                },
                c1: Fp2 {
                    c0: Fp(0x214204f9659218ef),
                    c1: Fp(0x96f25fd3ef81896),
                },
                c2: Fp2 {
                    c0: Fp(0x9251c878008fc0a),
                    c1: Fp(0x3e774da4e79e8a48),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x36832f29a238fb5a),
                    c1: Fp(0x1cf47a6fbba8b3ec),
                },
                c1: Fp2 {
                    c0: Fp(0x1c8953b3c44f6dfe),
                    c1: Fp(0x3bc825d1a69400e2),
                },
                c2: Fp2 {
                    c0: Fp(0x1e74c2cf8e10d7c),
                    c1: Fp(0xc49ef45b1e4117b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x38daa21a7587681b),
                    c1: Fp(0x233fd9db55b74389),
                },
                c1: Fp2 {
                    c0: Fp(0x21b317b399d3c6ff),
                    c1: Fp(0x188bbbf612f5b3bd),
                },
                c2: Fp2 {
                    c0: Fp(0x32d173b2a192cb04),
                    c1: Fp(0x131ef586933d2d2d),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3677011692915231),
                    c1: Fp(0x96916cbf186f19b),
                },
                c1: Fp2 {
                    c0: Fp(0x11a4c44efd12cdb1),
                    c1: Fp(0xb8ce70f2d41fc1e),
                },
                c2: Fp2 {
                    c0: Fp(0x33da11397a909fb6),
                    c1: Fp(0x1957d95f324efdee),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x339f1b9504a3c096),
                    c1: Fp(0x25c2a180eed5af82),
                },
                c1: Fp2 {
                    c0: Fp(0x8f628bd2110774b),
                    c1: Fp(0x31e1d9d4408c4e2b),
                },
                c2: Fp2 {
                    c0: Fp(0x167f245850f0e840),
                    c1: Fp(0x393b227116d1b751),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe51033821ceed3c),
                    c1: Fp(0xf5cb8bb9d9bea5a),
                },
                c1: Fp2 {
                    c0: Fp(0x3e0abe9f37f2d738),
                    c1: Fp(0x31a9f52cb21d9ff4),
                },
                c2: Fp2 {
                    c0: Fp(0x38bf29b3557ccf68),
                    c1: Fp(0x3bfc99414ae3621c),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x16c1113ca6e2672c),
                    c1: Fp(0x203e825672178dc6),
                },
                c1: Fp2 {
                    c0: Fp(0x31296e7499e6cf10),
                    c1: Fp(0x175ebe14e7be3dc7),
                },
                c2: Fp2 {
                    c0: Fp(0x1284d2d193f3df3d),
                    c1: Fp(0x1f723fe9e74e7d2d),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x13abfeb019cfc2b8),
                    c1: Fp(0x13a61667aa7477a8),
                },
                c1: Fp2 {
                    c0: Fp(0x33cbd707fd4017a),
                    c1: Fp(0xd9e7f1c5a1e8874),
                },
                c2: Fp2 {
                    c0: Fp(0x24611e1c9ce84ef7),
                    c1: Fp(0x1c3356a26dea5aab),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1c62e6c401683c98),
                    c1: Fp(0xfa3ab9afbfd10fa),
                },
                c1: Fp2 {
                    c0: Fp(0x257f0fb4679a1cdd),
                    c1: Fp(0x299a687cd6bed73a),
                },
                c2: Fp2 {
                    c0: Fp(0x8fc9135ca174863),
                    c1: Fp(0x2d9d47196c192df3),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2da229eb1531ffc),
                    c1: Fp(0x68dea68214ed2c3),
                },
                c1: Fp2 {
                    c0: Fp(0x3393aa941503ff05),
                    c1: Fp(0x3e2a527b21c0cfda),
                },
                c2: Fp2 {
                    c0: Fp(0x1d1e56243f39145f),
                    c1: Fp(0x111f6305f9b57400),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x51546cb5ef15e40),
                    c1: Fp(0x32519911bfed8aa),
                },
                c1: Fp2 {
                    c0: Fp(0x1764905297fe5db6),
                    c1: Fp(0x3e8548b1c1c985b3),
                },
                c2: Fp2 {
                    c0: Fp(0x5e287ddc204de8a),
                    c1: Fp(0x1e59ac469fc82b48),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2e36043e0c8865f3),
                    c1: Fp(0xb30eb1c5f63f05d),
                },
                c1: Fp2 {
                    c0: Fp(0x12b02a52517471de),
                    c1: Fp(0x41565e3d272992ff),
                },
                c2: Fp2 {
                    c0: Fp(0x3034d5101ec47c80),
                    c1: Fp(0x2046465660e73d63),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x8350c8d80c140c4),
                    c1: Fp(0x386868beecf6dade),
                },
                c1: Fp2 {
                    c0: Fp(0x1ad04df2817e1223),
                    c1: Fp(0xbdb4bf6c54b2481),
                },
                c2: Fp2 {
                    c0: Fp(0x323bdc7580f2b2a2),
                    c1: Fp(0x55b2ad9d360a4ca),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x28b9aace73a19876),
                    c1: Fp(0x1c0e701e46e1fcf),
                },
                c1: Fp2 {
                    c0: Fp(0xd6ed978f528d41c),
                    c1: Fp(0x1642cdc175b2b8c4),
                },
                c2: Fp2 {
                    c0: Fp(0xdfec68934ecb05b),
                    c1: Fp(0x14fbaaf990bf872),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a7e95413371f6a2),
                    c1: Fp(0x2e5ad901ed568466),
                },
                c1: Fp2 {
                    c0: Fp(0x80ba2c271182758),
                    c1: Fp(0x1a69b2e439bd749e),
                },
                c2: Fp2 {
                    c0: Fp(0x288ce70a26c7676b),
                    c1: Fp(0x116d8976a8bf7f13),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2b5d421c05b221f9),
                    c1: Fp(0x4028bbcdcb23b78c),
                },
                c1: Fp2 {
                    c0: Fp(0x308b31b2c5c76615),
                    c1: Fp(0x2f35b9dbaef8a049),
                },
                c2: Fp2 {
                    c0: Fp(0x1b49b700251f60ff),
                    c1: Fp(0x407997dbf7cbfb8d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2c166e399a09ed57),
                    c1: Fp(0x359ef79094c23464),
                },
                c1: Fp2 {
                    c0: Fp(0x34b3c749b11d8e2c),
                    c1: Fp(0x14e2ffd119b5926d),
                },
                c2: Fp2 {
                    c0: Fp(0x203d7327b9d4ef94),
                    c1: Fp(0x192b22684dd9d7a9),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x32628ad2520d21e2),
                    c1: Fp(0xcc340991332aca2),
                },
                c1: Fp2 {
                    c0: Fp(0x68c927e1fbd4558),
                    c1: Fp(0xc91acd44c4c6a5b),
                },
                c2: Fp2 {
                    c0: Fp(0x3d740c9d1d81386d),
                    c1: Fp(0x25ad5d3786a44fbd),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x36cd1312aba6e9d8),
                    c1: Fp(0x3940b9f2cffb30d7),
                },
                c1: Fp2 {
                    c0: Fp(0x162a7baa6b8201d2),
                    c1: Fp(0x137a2b4b66b10432),
                },
                c2: Fp2 {
                    c0: Fp(0x2d865e1a3883bb8a),
                    c1: Fp(0x3f743d99e0399de8),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2489dd206829a7b4),
                    c1: Fp(0x2c8a14de35c107bd),
                },
                c1: Fp2 {
                    c0: Fp(0x1f49670ca4b4c859),
                    c1: Fp(0x1492838d13316389),
                },
                c2: Fp2 {
                    c0: Fp(0x3c650afd98b018f4),
                    c1: Fp(0x257880f95b28b14b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x194e8dd5cdf03035),
                    c1: Fp(0x3bd547921b39ee18),
                },
                c1: Fp2 {
                    c0: Fp(0x2d5d2cc1c38d534d),
                    c1: Fp(0x2e33055381e74066),
                },
                c2: Fp2 {
                    c0: Fp(0x3287e13910fb5e1f),
                    c1: Fp(0x22bbe8d66701185),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xe8a92a04f649864),
                    c1: Fp(0x3387bf4645b7886a),
                },
                c1: Fp2 {
                    c0: Fp(0x1872dc347bbb6874),
                    c1: Fp(0x1ba3c9d782dc1597),
                },
                c2: Fp2 {
                    c0: Fp(0x14b6203890c70d99),
                    c1: Fp(0x40049459be6de44),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a42fc12410622b4),
                    c1: Fp(0x2d17d6e1e5dbf9aa),
                },
                c1: Fp2 {
                    c0: Fp(0x2bd6aa7bb2dc9e9e),
                    c1: Fp(0x1ea22fb82d97978),
                },
                c2: Fp2 {
                    c0: Fp(0x11e4571166a7823f),
                    c1: Fp(0x220b6ff78602785b),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x23ff99aa62c3039d),
                    c1: Fp(0x149446580668c6fa),
                },
                c1: Fp2 {
                    c0: Fp(0x313170242799d522),
                    c1: Fp(0x26839707612e7d78),
                },
                c2: Fp2 {
                    c0: Fp(0x1701bb630cb1ba3f),
                    c1: Fp(0xdf00cf7a0baffce),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3a7545ccb71ece79),
                    c1: Fp(0x388f44d9c9789169),
                },
                c1: Fp2 {
                    c0: Fp(0x30ded7e070e70ead),
                    c1: Fp(0x122e3299fe80ab2a),
                },
                c2: Fp2 {
                    c0: Fp(0x113f437b0bd25a1b),
                    c1: Fp(0x9906df7d927e27f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x22484bcedce7b519),
                    c1: Fp(0x3f479875c38e0a9c),
                },
                c1: Fp2 {
                    c0: Fp(0x202c1c941a4da81f),
                    c1: Fp(0x2c384369de2d42ac),
                },
                c2: Fp2 {
                    c0: Fp(0x8dfb23054a6ad49),
                    c1: Fp(0x258828e914ca4ef9),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xc3697e2daf9a0),
                    c1: Fp(0x98ceff4be15e431),
                },
                c1: Fp2 {
                    c0: Fp(0x190f8591aa09493a),
                    c1: Fp(0x2efed4f9b89ddd8a),
                },
                c2: Fp2 {
                    c0: Fp(0xe6214f225746f79),
                    c1: Fp(0x398ff060b09e4460),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x238156baf44da087),
                    c1: Fp(0x24553aee7973d0a3),
                },
                c1: Fp2 {
                    c0: Fp(0x1c8d3ed721d4676f),
                    c1: Fp(0xcfe9e1d42a2ea61),
                },
                c2: Fp2 {
                    c0: Fp(0x40ee4e109d402286),
                    c1: Fp(0x3c244b03b3c79519),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2903d2367bd22404),
                    c1: Fp(0x39bfcb650ae46484),
                },
                c1: Fp2 {
                    c0: Fp(0x1e4b9fbee65444ed),
                    c1: Fp(0x2f9e3c36b4d863bc),
                },
                c2: Fp2 {
                    c0: Fp(0x3669d52c9b89cd9f),
                    c1: Fp(0x1c9ece4bd0baa142),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2e35473be2cb86c4),
                    c1: Fp(0xd322ebc3cffca6c),
                },
                c1: Fp2 {
                    c0: Fp(0x4cdbddea9ee2f6f),
                    c1: Fp(0x21e9ea572a61b0e5),
                },
                c2: Fp2 {
                    c0: Fp(0x10737a464ea3497d),
                    c1: Fp(0x29b0c351f6282511),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2e57e6d746d58229),
                    c1: Fp(0x33ce12a56497bf72),
                },
                c1: Fp2 {
                    c0: Fp(0x335ec420e1413f82),
                    c1: Fp(0x14f565f5ccc1bcf7),
                },
                c2: Fp2 {
                    c0: Fp(0x2df680893fd5218f),
                    c1: Fp(0x398799b490e1b44),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1a4bf8a78c8ba5a8),
                    c1: Fp(0x314423c5ada98161),
                },
                c1: Fp2 {
                    c0: Fp(0x2237bcdb08b1f2d7),
                    c1: Fp(0x86e78b6a5c7fb10),
                },
                c2: Fp2 {
                    c0: Fp(0x35fd796c45ed0ace),
                    c1: Fp(0x34f7fea25070b77e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x18877d242ae0713d),
                    c1: Fp(0x1933c72d1eef8478),
                },
                c1: Fp2 {
                    c0: Fp(0x3cc6aa12a628a988),
                    c1: Fp(0x1c4fbf0e28614dce),
                },
                c2: Fp2 {
                    c0: Fp(0x15688dc178ed3951),
                    c1: Fp(0x1c52695f788675b9),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x30c80420ed38bd68),
                    c1: Fp(0x976207365420),
                },
                c1: Fp2 {
                    c0: Fp(0x2ada9bcc10a368e),
                    c1: Fp(0x3aafc304506b0287),
                },
                c2: Fp2 {
                    c0: Fp(0x1b9792d80a3e4d81),
                    c1: Fp(0x1e19592dcf7f6575),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1b71e78e243876a5),
                    c1: Fp(0x34049bed26cf7880),
                },
                c1: Fp2 {
                    c0: Fp(0x395ece3b5acbc5dc),
                    c1: Fp(0x1d495ea0da793a49),
                },
                c2: Fp2 {
                    c0: Fp(0x40d21f51086ce6d4),
                    c1: Fp(0x335b82cabf0df267),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x20acf0c4371c65),
                    c1: Fp(0x30a7c1240853dda2),
                },
                c1: Fp2 {
                    c0: Fp(0x3162f983d18f6165),
                    c1: Fp(0x1102d9c223b1742c),
                },
                c2: Fp2 {
                    c0: Fp(0x40f4ca21901fe96),
                    c1: Fp(0x581c18e247077fc),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3db2899d0137ad76),
                    c1: Fp(0x21d7be085780318a),
                },
                c1: Fp2 {
                    c0: Fp(0x67890888ee03f69),
                    c1: Fp(0x13286a35919b1d21),
                },
                c2: Fp2 {
                    c0: Fp(0x12aad59636103292),
                    c1: Fp(0x1246cd09f9d3763f),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1885fedc7bffef08),
                    c1: Fp(0x32eaa2240f1a658e),
                },
                c1: Fp2 {
                    c0: Fp(0xd16c8d2bf78c53f),
                    c1: Fp(0x2e76cd0d049c6a1e),
                },
                c2: Fp2 {
                    c0: Fp(0x1aca92462104a729),
                    c1: Fp(0x226bffb3d6407f49),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1a49b020f6aec0e9),
                    c1: Fp(0x3114e551bc70c500),
                },
                c1: Fp2 {
                    c0: Fp(0x75c4fbe4cd84984),
                    c1: Fp(0x553a4f4a63f8faa),
                },
                c2: Fp2 {
                    c0: Fp(0x129116ef735474a5),
                    c1: Fp(0x233ae8c47bdbd909),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3f105b86bf3eee2),
                    c1: Fp(0x304253ecc36ebb06),
                },
                c1: Fp2 {
                    c0: Fp(0x122982c8028e29a6),
                    c1: Fp(0x336d21386fbb2b91),
                },
                c2: Fp2 {
                    c0: Fp(0x83767702d44bee9),
                    c1: Fp(0x1c435cea64183cf0),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x30831459ef8ba879),
                    c1: Fp(0x2c4a617a69353a81),
                },
                c1: Fp2 {
                    c0: Fp(0x268d64bd2e712816),
                    c1: Fp(0x1eb6c31eccd7746a),
                },
                c2: Fp2 {
                    c0: Fp(0x109dd7eafd91d44d),
                    c1: Fp(0x15f4bce44a1e694a),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x298c5703ab8a9f87),
                    c1: Fp(0x16560f6acb8a1c1a),
                },
                c1: Fp2 {
                    c0: Fp(0x251cd20583e571fe),
                    c1: Fp(0x3a3ba81d49f19a09),
                },
                c2: Fp2 {
                    c0: Fp(0x2e8f6c73d034188f),
                    c1: Fp(0x2da0307c1cb2c6e4),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2fbdff6180f842d0),
                    c1: Fp(0xcfc686dd0933886),
                },
                c1: Fp2 {
                    c0: Fp(0x3dc362a29a2ade8e),
                    c1: Fp(0x1be3864fadad7e43),
                },
                c2: Fp2 {
                    c0: Fp(0x2e268afa3d99528),
                    c1: Fp(0x2af119f67076d49e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x19e04037cf0b9666),
                    c1: Fp(0x2cc59fa060da083d),
                },
                c1: Fp2 {
                    c0: Fp(0x1e2e40505dc0b454),
                    c1: Fp(0x3aaf7d80416184e8),
                },
                c2: Fp2 {
                    c0: Fp(0x1d6acfcc7e10c070),
                    c1: Fp(0x1dad80c124ea6335),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3ff0b89e7e5f09d8),
                    c1: Fp(0x3ce0a1275b9ef997),
                },
                c1: Fp2 {
                    c0: Fp(0x47c9aac3c1f0b73),
                    c1: Fp(0x626d9cd3057a5bb),
                },
                c2: Fp2 {
                    c0: Fp(0x2999217c5647753b),
                    c1: Fp(0xd03d54e2fe029c6),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x26ee2f42252a861b),
                    c1: Fp(0x31ff84030460ab4b),
                },
                c1: Fp2 {
                    c0: Fp(0x1dbc092d8e9b7e78),
                    c1: Fp(0x243f7b1229aef2b7),
                },
                c2: Fp2 {
                    c0: Fp(0x1d9255cecd88dac5),
                    c1: Fp(0x2d149b327cb94f4e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xbb3ac3d631a967f),
                    c1: Fp(0x1c1492a3b1cd63b7),
                },
                c1: Fp2 {
                    c0: Fp(0xcdb961683e56c29),
                    c1: Fp(0x108559e255b61ee0),
                },
                c2: Fp2 {
                    c0: Fp(0xd49efbb3da11c2f),
                    c1: Fp(0x32d7bcb442f27ca6),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x13960ae8acb9334a),
                    c1: Fp(0x333f5e2755a35add),
                },
                c1: Fp2 {
                    c0: Fp(0x27f4a8cf98fb84e8),
                    c1: Fp(0x3df5e2f32aee5ee8),
                },
                c2: Fp2 {
                    c0: Fp(0x129655ed1e90626),
                    c1: Fp(0x2feab5f707c5f5e3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xb7168ed23158239),
                    c1: Fp(0x2b3b7a0e5f825e6b),
                },
                c1: Fp2 {
                    c0: Fp(0x2ba1f7a58ff55045),
                    c1: Fp(0x3be18f5eac73140d),
                },
                c2: Fp2 {
                    c0: Fp(0x97f946cdad97614),
                    c1: Fp(0x32ad36e5a11dd0bb),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3f228af76de0a513),
                    c1: Fp(0xbb64060057e42d0),
                },
                c1: Fp2 {
                    c0: Fp(0x25ccfa229b59cab5),
                    c1: Fp(0x2baa6753584e71f4),
                },
                c2: Fp2 {
                    c0: Fp(0x96b456e97b603ed),
                    c1: Fp(0x1d3e5f5a37b43e4f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x264da821f09f4b98),
                    c1: Fp(0x379f68d586418635),
                },
                c1: Fp2 {
                    c0: Fp(0x12271df129e85717),
                    c1: Fp(0x151cef4f5e69976f),
                },
                c2: Fp2 {
                    c0: Fp(0xa551999e4546db6),
                    c1: Fp(0x233a1a5e0317dd9),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xba7754bcf844a62),
                    c1: Fp(0x1e6b0e26ded14bf3),
                },
                c1: Fp2 {
                    c0: Fp(0x27a4ec18d8e76d20),
                    c1: Fp(0x26bf362b4229e5c3),
                },
                c2: Fp2 {
                    c0: Fp(0x2cc7163bb391c34d),
                    c1: Fp(0x195057115a87bd76),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x101bb9abccdcd004),
                    c1: Fp(0x8c28b1d3b2e3589),
                },
                c1: Fp2 {
                    c0: Fp(0x289ccb308fcc9965),
                    c1: Fp(0x26c4cb0ff513c4b),
                },
                c2: Fp2 {
                    c0: Fp(0x1eb057c53454f16d),
                    c1: Fp(0x24c61468f20d7171),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x11fc734a99b54d33),
                    c1: Fp(0x38a02baa75117ed7),
                },
                c1: Fp2 {
                    c0: Fp(0x10ee540a68bde3a1),
                    c1: Fp(0x3059b1d666b44276),
                },
                c2: Fp2 {
                    c0: Fp(0x1c0527195bcdebb9),
                    c1: Fp(0x7d5e74b48984783),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1550d8a556659284),
                    c1: Fp(0x14ea156e332e1eb4),
                },
                c1: Fp2 {
                    c0: Fp(0x278eee2176c8f03f),
                    c1: Fp(0x3ceb3ba5d7a9370),
                },
                c2: Fp2 {
                    c0: Fp(0x2fab99053eaee7b5),
                    c1: Fp(0x2ad0124462196849),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c447161a58948d2),
                    c1: Fp(0x195c848702a2e961),
                },
                c1: Fp2 {
                    c0: Fp(0x2d0fefd8708ef19e),
                    c1: Fp(0x4def066079d77ff),
                },
                c2: Fp2 {
                    c0: Fp(0x1af628601b354e2d),
                    c1: Fp(0xdea50562bc1b5f8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x417437cdaa01a3e7),
                    c1: Fp(0x25f1a1202aa72f9e),
                },
                c1: Fp2 {
                    c0: Fp(0x2286737d479c66bb),
                    c1: Fp(0x3fb837875d69a6e6),
                },
                c2: Fp2 {
                    c0: Fp(0x243dbb3cf90c81c3),
                    c1: Fp(0x37edad07ce1b6759),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x368247298834ca83),
                    c1: Fp(0x21546606ee7266e1),
                },
                c1: Fp2 {
                    c0: Fp(0x1fc5151c3b2d0f83),
                    c1: Fp(0x6cfccdaec0bc917),
                },
                c2: Fp2 {
                    c0: Fp(0x838efc353227eba),
                    c1: Fp(0x228ec50ea142351a),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x18b714117cc02c49),
                    c1: Fp(0xdd83da97c311eb3),
                },
                c1: Fp2 {
                    c0: Fp(0x19f45df898ca6736),
                    c1: Fp(0x10fa64f90ff0e50e),
                },
                c2: Fp2 {
                    c0: Fp(0x199ea289f13aae1c),
                    c1: Fp(0x2cbdc46d2dc32ac2),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x32740a385325379f),
                    c1: Fp(0x7d1cee384bdfc8b),
                },
                c1: Fp2 {
                    c0: Fp(0x393452c7a162d332),
                    c1: Fp(0x19197b4e146496b4),
                },
                c2: Fp2 {
                    c0: Fp(0x268ca00f243e9bc8),
                    c1: Fp(0x1ccccb399b367adf),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x333e2c9afcd99afd),
                    c1: Fp(0x3854cee816d0d589),
                },
                c1: Fp2 {
                    c0: Fp(0x3e4ba95cf96cfc00),
                    c1: Fp(0x2912d5338596f470),
                },
                c2: Fp2 {
                    c0: Fp(0x38b04e44614af18b),
                    c1: Fp(0x2df7eb0eb7b0eb59),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xbed64bf863b057d),
                    c1: Fp(0x11e5e79cb644f9b9),
                },
                c1: Fp2 {
                    c0: Fp(0x3cce31203736251c),
                    c1: Fp(0xe45b0c6a072209),
                },
                c2: Fp2 {
                    c0: Fp(0x1027af59a6dffbaa),
                    c1: Fp(0x109b7ebfa43c8e40),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x5600bf0abdd68d6),
                    c1: Fp(0x3904d80684cb3844),
                },
                c1: Fp2 {
                    c0: Fp(0x37285a8ebc134c58),
                    c1: Fp(0x8d33f63e1909f69),
                },
                c2: Fp2 {
                    c0: Fp(0xaf3fcc79ac9ca0),
                    c1: Fp(0x3553866b01f909c1),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x40fdc8766ed64913),
                    c1: Fp(0x7e19c374467a05b),
                },
                c1: Fp2 {
                    c0: Fp(0x19a716180b302bbe),
                    c1: Fp(0xf5ed833f424fe3e),
                },
                c2: Fp2 {
                    c0: Fp(0x229cd15c4edbbe34),
                    c1: Fp(0x69ae61b8621acaf),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1198f150d6c8e115),
                    c1: Fp(0x8e02758de59fee8),
                },
                c1: Fp2 {
                    c0: Fp(0x20e5e224e9bc4707),
                    c1: Fp(0x2bd749c5c31ec123),
                },
                c2: Fp2 {
                    c0: Fp(0x1174b0697cf14600),
                    c1: Fp(0x3020d10837b80ba4),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x7b42ebe55dd860d),
                    c1: Fp(0x397ebba150cde809),
                },
                c1: Fp2 {
                    c0: Fp(0x3e053f238fa99e2b),
                    c1: Fp(0x2b1766a1b5e2f0d3),
                },
                c2: Fp2 {
                    c0: Fp(0x24a4e53dafe77f93),
                    c1: Fp(0xbd85a8fb9a89b8),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd9e87df488244ee),
                    c1: Fp(0x1b7182e5f30730a2),
                },
                c1: Fp2 {
                    c0: Fp(0x908bb6787dd77aa),
                    c1: Fp(0x2c7901f8e4e78460),
                },
                c2: Fp2 {
                    c0: Fp(0xbb44f3e12232c16),
                    c1: Fp(0x3f00784338663bca),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xc209cc5eacd9c40),
                    c1: Fp(0x1640247af81184fa),
                },
                c1: Fp2 {
                    c0: Fp(0x3271de38fd39f48e),
                    c1: Fp(0x150b67ea6f95af3),
                },
                c2: Fp2 {
                    c0: Fp(0x37ee53e144eb7bae),
                    c1: Fp(0x3c91ab96a52bcbc2),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x249ea74c6ae0ee49),
                    c1: Fp(0x39670d73804c4792),
                },
                c1: Fp2 {
                    c0: Fp(0x128fb1fff24e4fd0),
                    c1: Fp(0x13ba7de25858b847),
                },
                c2: Fp2 {
                    c0: Fp(0x3ac019d112cb2661),
                    c1: Fp(0x18f19c9cc9f696b7),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x9550b78148d2b7d),
                    c1: Fp(0x9b436691e34c365),
                },
                c1: Fp2 {
                    c0: Fp(0x78a5e9736bdd4a6),
                    c1: Fp(0x1fb7e105edf4767e),
                },
                c2: Fp2 {
                    c0: Fp(0xb805f859e959b99),
                    c1: Fp(0x49d5c1da4cdbd0e),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x156b68542c3308b8),
                    c1: Fp(0x1936aa3864ab27eb),
                },
                c1: Fp2 {
                    c0: Fp(0x2c3ee06eca4a3d36),
                    c1: Fp(0x298fee87ffbe0625),
                },
                c2: Fp2 {
                    c0: Fp(0x18e566ad52ea8ef7),
                    c1: Fp(0x8f60380304e721b),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x606aaff14dfafe7),
                    c1: Fp(0x3304f76beabfb246),
                },
                c1: Fp2 {
                    c0: Fp(0xed1356ec8e8e23e),
                    c1: Fp(0x3f4c39996be798b4),
                },
                c2: Fp2 {
                    c0: Fp(0x12499b566700d6f7),
                    c1: Fp(0x1db3e00e06a91e87),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd1816fa651d0c2c),
                    c1: Fp(0x16569073397501b2),
                },
                c1: Fp2 {
                    c0: Fp(0x189da861c4865826),
                    c1: Fp(0x1edd5205f7e79721),
                },
                c2: Fp2 {
                    c0: Fp(0x1291e27fc8504618),
                    c1: Fp(0x238c88a4c37fffe1),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2dd45a26a71cee3f),
                    c1: Fp(0x304bd46c8e186fe3),
                },
                c1: Fp2 {
                    c0: Fp(0x29716a237a697bdc),
                    c1: Fp(0x1152947da76d982f),
                },
                c2: Fp2 {
                    c0: Fp(0x11caabfc5fa441d5),
                    c1: Fp(0x28a06860c76d0eef),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x22abe82d24c0adf5),
                    c1: Fp(0x2ad1f5d2b9d48688),
                },
                c1: Fp2 {
                    c0: Fp(0x57463dec89b80ac),
                    c1: Fp(0x1f544da5ac99e574),
                },
                c2: Fp2 {
                    c0: Fp(0x24fe87b9c7cfd9f5),
                    c1: Fp(0x56ce7e9b57c7325),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x35bdaca272a06f63),
                    c1: Fp(0x2c96f8a83ce27341),
                },
                c1: Fp2 {
                    c0: Fp(0x4068d74b6b0d4030),
                    c1: Fp(0x81b50081d0a3d8d),
                },
                c2: Fp2 {
                    c0: Fp(0x33c99a1c38e164d9),
                    c1: Fp(0x212d97c56b708169),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x388c22c66751e36),
                    c1: Fp(0x3182f9dc566b8f3f),
                },
                c1: Fp2 {
                    c0: Fp(0x1540b2744d0010b6),
                    c1: Fp(0x39948fc48192ddd2),
                },
                c2: Fp2 {
                    c0: Fp(0x3d01da1b46578cd9),
                    c1: Fp(0x109e51fe4b202011),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1f3ab21d02c977da),
                    c1: Fp(0x17d82e718af69b56),
                },
                c1: Fp2 {
                    c0: Fp(0x6a0414382a56d48),
                    c1: Fp(0x14319ff9cb9df5da),
                },
                c2: Fp2 {
                    c0: Fp(0x391bee81814f379c),
                    c1: Fp(0x353c5e2c3ed417c6),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xc2cb6a95a0b872e),
                    c1: Fp(0x2ecfaef49f929f9f),
                },
                c1: Fp2 {
                    c0: Fp(0x190b0210f7e12ce5),
                    c1: Fp(0x228765b8e943dd17),
                },
                c2: Fp2 {
                    c0: Fp(0x1d04b8ff3380d597),
                    c1: Fp(0x31f40a92e96a6d76),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3720cb8299a79b5),
                    c1: Fp(0x639492b858cf7ae),
                },
                c1: Fp2 {
                    c0: Fp(0x33a1700ac24d5a49),
                    c1: Fp(0xf67387320869d27),
                },
                c2: Fp2 {
                    c0: Fp(0x9cad2ea76506087),
                    c1: Fp(0xb6852801cf6b080),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x37ccfbb02c9de3ce),
                    c1: Fp(0x2d7309d52edd1e3a),
                },
                c1: Fp2 {
                    c0: Fp(0x185ce204ca1274f),
                    c1: Fp(0x34c480ef8fc3bc11),
                },
                c2: Fp2 {
                    c0: Fp(0x1fcaac945ae7f48b),
                    c1: Fp(0x30f9a8d0831b7001),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x15c97af294a63a4c),
                    c1: Fp(0x3557146e3338686e),
                },
                c1: Fp2 {
                    c0: Fp(0x2f0e4fe69a74df10),
                    c1: Fp(0x27b4053f6da06ad2),
                },
                c2: Fp2 {
                    c0: Fp(0x2af457cf671225dd),
                    c1: Fp(0x38a86ffb5160a49b),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3d41524849967b0d),
                    c1: Fp(0x16481c2e703832c4),
                },
                c1: Fp2 {
                    c0: Fp(0x9e60ca74b8b009e),
                    c1: Fp(0x379b5097c9b65bfb),
                },
                c2: Fp2 {
                    c0: Fp(0x3a3e16b9504918ac),
                    c1: Fp(0x1dd02216791e89be),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x32a71847c6b4c6dd),
                    c1: Fp(0x2502f23d90e6d2a2),
                },
                c1: Fp2 {
                    c0: Fp(0x130ce246def8a737),
                    c1: Fp(0x27bf3f65ffdb86b8),
                },
                c2: Fp2 {
                    c0: Fp(0x2f4b2fcd5d76310b),
                    c1: Fp(0x3e0b2f2863a7a7a7),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x18183394906c0555),
                    c1: Fp(0x301cf2371684be20),
                },
                c1: Fp2 {
                    c0: Fp(0x1ddcb56db48c02f1),
                    c1: Fp(0x124a5d08648a4adf),
                },
                c2: Fp2 {
                    c0: Fp(0x3a78e347c57f6778),
                    c1: Fp(0x3fa2b6fa657ca851),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xbd9c8769b189592),
                    c1: Fp(0x25beddeada70688d),
                },
                c1: Fp2 {
                    c0: Fp(0x382ed22dfbeec09),
                    c1: Fp(0x6b5d46122f5046d),
                },
                c2: Fp2 {
                    c0: Fp(0x25e5365279ffda62),
                    c1: Fp(0x2fcfa71ce00a5832),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3d3484e87a21f086),
                    c1: Fp(0x5b398bd4ee2d101),
                },
                c1: Fp2 {
                    c0: Fp(0x1ffc2f1d5d25667d),
                    c1: Fp(0x1b0f959b16a4f00f),
                },
                c2: Fp2 {
                    c0: Fp(0x100455b6a1c45899),
                    c1: Fp(0x40fdad565a93bb33),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x170cbec0c05f954),
                    c1: Fp(0x3c16ce90b0288209),
                },
                c1: Fp2 {
                    c0: Fp(0x2d35882016640f3e),
                    c1: Fp(0x3f51e54630147e45),
                },
                c2: Fp2 {
                    c0: Fp(0xdd1b181d8510de7),
                    c1: Fp(0x281c64e73a8c061d),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x34f42fc4ff7633ab),
                    c1: Fp(0x7d68ebfb74584f7),
                },
                c1: Fp2 {
                    c0: Fp(0x37179a996d4fbbb3),
                    c1: Fp(0x1f109ab586a9835b),
                },
                c2: Fp2 {
                    c0: Fp(0x189ebb626a8de525),
                    c1: Fp(0x54296c588e3ebc3),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c9a0b705aac3b32),
                    c1: Fp(0x26ea3c153077a8a0),
                },
                c1: Fp2 {
                    c0: Fp(0x2d22eeb5f273f6e5),
                    c1: Fp(0x4090dad5bd16170d),
                },
                c2: Fp2 {
                    c0: Fp(0x3683a1bdaa5ae28c),
                    c1: Fp(0x3e59bdb568181c21),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x16aff1d9427245c9),
                    c1: Fp(0x16f861748365a0d4),
                },
                c1: Fp2 {
                    c0: Fp(0x2e3b639650a7c6d1),
                    c1: Fp(0x3aa1fdba83344ebd),
                },
                c2: Fp2 {
                    c0: Fp(0x5eb1a44a67cbf69),
                    c1: Fp(0x35255a18043aa4d8),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2473a881dd3e715e),
                    c1: Fp(0x112d2b4b3b1ed218),
                },
                c1: Fp2 {
                    c0: Fp(0x147439d42be1f4b5),
                    c1: Fp(0x6017aed531f074b),
                },
                c2: Fp2 {
                    c0: Fp(0xcd0f770105ea5db),
                    c1: Fp(0x1dda85869bc23cd),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x5ee1cb206275af),
                    c1: Fp(0x3dcbe9ecd89cb8af),
                },
                c1: Fp2 {
                    c0: Fp(0x1ce71ec2bcffd461),
                    c1: Fp(0x1e9a61a6c1bca270),
                },
                c2: Fp2 {
                    c0: Fp(0xd91c18829d5b773),
                    c1: Fp(0x272931f47055a987),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x148b2a297dd2f779),
                    c1: Fp(0x3ca9bc8bd49a1f4a),
                },
                c1: Fp2 {
                    c0: Fp(0x15046e258374ba96),
                    c1: Fp(0x1737275f277166da),
                },
                c2: Fp2 {
                    c0: Fp(0x357e1ba19c000cdc),
                    c1: Fp(0x31cca3c6a48ecdca),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3e698eceb147d48e),
                    c1: Fp(0x2c9527b8e57b0e61),
                },
                c1: Fp2 {
                    c0: Fp(0x5bf40328c3521b9),
                    c1: Fp(0x400181cca60ebc6e),
                },
                c2: Fp2 {
                    c0: Fp(0x8638365706f5d2e),
                    c1: Fp(0x11fea0366f823b12),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x33f094b83b6b02d4),
                    c1: Fp(0x3c9f2dd77fa633d2),
                },
                c1: Fp2 {
                    c0: Fp(0x4b531e468207d8a),
                    c1: Fp(0x693d07b374dc77b),
                },
                c2: Fp2 {
                    c0: Fp(0xbe4c71ed6556b1e),
                    c1: Fp(0x31f152dc203f47ca),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x11833253b0586e71),
                    c1: Fp(0x3e9aaa1c9de47408),
                },
                c1: Fp2 {
                    c0: Fp(0x3892bfba198d5085),
                    c1: Fp(0x8b47b2ecfcac2d),
                },
                c2: Fp2 {
                    c0: Fp(0xdfb1881ba3e89bf),
                    c1: Fp(0x14a0bec3d241f9b7),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ca4d1600a681008),
                    c1: Fp(0x2aef230f21569475),
                },
                c1: Fp2 {
                    c0: Fp(0x3d26a87408034867),
                    c1: Fp(0x75a36c9c4091883),
                },
                c2: Fp2 {
                    c0: Fp(0x18b9257c578fd8cf),
                    c1: Fp(0x7d170bde2075661),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3487dc71c88f817e),
                    c1: Fp(0x21be840964ed8c47),
                },
                c1: Fp2 {
                    c0: Fp(0x1cde76fbb1fbfff0),
                    c1: Fp(0x26144a8ce4d83aa0),
                },
                c2: Fp2 {
                    c0: Fp(0x98c4b81a3d67d63),
                    c1: Fp(0x286b1f67624bb3d0),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3c4233d12ad4b862),
                    c1: Fp(0xa86709d8d4c20b9),
                },
                c1: Fp2 {
                    c0: Fp(0x3b9b1c03962800d2),
                    c1: Fp(0xcad60d69703e4e2),
                },
                c2: Fp2 {
                    c0: Fp(0x2cd7275e2643eb1b),
                    c1: Fp(0x14696d0aa96b7a22),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xa3dd1b66f41712f),
                    c1: Fp(0x2129cb67a352a831),
                },
                c1: Fp2 {
                    c0: Fp(0xa98b41d7dc1cfcd),
                    c1: Fp(0x1bbc5323c95fae0d),
                },
                c2: Fp2 {
                    c0: Fp(0x1f3bb04564144aa1),
                    c1: Fp(0x12fe5be0f5545696),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x25fe20099c4d64e4),
                    c1: Fp(0x16ecc72be08a3fbf),
                },
                c1: Fp2 {
                    c0: Fp(0x6316a71c7ea7280),
                    c1: Fp(0x2fdff316bb40ad25),
                },
                c2: Fp2 {
                    c0: Fp(0x3cf8e2a9829a045c),
                    c1: Fp(0x34f6328e629cbca),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x9d48aabf1d295df),
                    c1: Fp(0x13a32bf78bc0cee9),
                },
                c1: Fp2 {
                    c0: Fp(0x1ec4dd011fe78dba),
                    c1: Fp(0x28b5d5009b8a9e18),
                },
                c2: Fp2 {
                    c0: Fp(0x32822e686bbffc04),
                    c1: Fp(0x329c266a04e50fb8),
                },
            },
            z: Fp6::one(),
        },
    ]),
    LookupTable([
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3d1bc5730340dc08),
                    c1: Fp(0x345cdb13a36e8c19),
                },
                c1: Fp2 {
                    c0: Fp(0x1a9a90fd64e072e2),
                    c1: Fp(0x6e1800d697504f2),
                },
                c2: Fp2 {
                    c0: Fp(0x24149077434d9ffc),
                    c1: Fp(0x8975955a5d4d2e3),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3aec8ccf9bcf9a8),
                    c1: Fp(0x2f243c7d873356b4),
                },
                c1: Fp2 {
                    c0: Fp(0x22839f94d18944c8),
                    c1: Fp(0x231913fa6f1473d4),
                },
                c2: Fp2 {
                    c0: Fp(0x36e8e1899de06f9e),
                    c1: Fp(0x3402c57fd4ed9af5),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1ff517172dc57903),
                    c1: Fp(0x1bd760ddeca99e8d),
                },
                c1: Fp2 {
                    c0: Fp(0x1eae008f236245d2),
                    c1: Fp(0x30b8e16a6616e511),
                },
                c2: Fp2 {
                    c0: Fp(0x20ac1e67446fff14),
                    c1: Fp(0x3609dd2c9971a447),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2580097296229300),
                    c1: Fp(0x4014e8674e82a0d8),
                },
                c1: Fp2 {
                    c0: Fp(0x1d19b14fca37a5e7),
                    c1: Fp(0x1b6155ed64028aa1),
                },
                c2: Fp2 {
                    c0: Fp(0x49a31a16a7fb1aa),
                    c1: Fp(0x3f1e5793258b467c),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x348f66c17fab6a77),
                    c1: Fp(0x32d67db2cf1ec556),
                },
                c1: Fp2 {
                    c0: Fp(0x220776a831140eeb),
                    c1: Fp(0x10425ed8f70a948a),
                },
                c2: Fp2 {
                    c0: Fp(0x22dc198d70245820),
                    c1: Fp(0x19725851db8b013f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x331af0b2bc655bd3),
                    c1: Fp(0x4324be4b1f25c5d),
                },
                c1: Fp2 {
                    c0: Fp(0x351dfcc064a76797),
                    c1: Fp(0x1a1030e098321178),
                },
                c2: Fp2 {
                    c0: Fp(0x3c073c980d401646),
                    c1: Fp(0x3bb7eb6da577358e),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x35a3edd56105636e),
                    c1: Fp(0x6267cb69c5fbfa9),
                },
                c1: Fp2 {
                    c0: Fp(0x396d494b87c3082f),
                    c1: Fp(0x1f3942194d0a9c38),
                },
                c2: Fp2 {
                    c0: Fp(0x101f0b2f6e7905e2),
                    c1: Fp(0x2c65702605ce7fe1),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x27dca52de9d41a41),
                    c1: Fp(0xdd770664a374f6e),
                },
                c1: Fp2 {
                    c0: Fp(0x193049d9f0936070),
                    c1: Fp(0x2d19dce54698eaf7),
                },
                c2: Fp2 {
                    c0: Fp(0x1c71f982651e59ae),
                    c1: Fp(0x2b4d37cb11f134e8),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x1e1f6d981de69c37),
                    c1: Fp(0x1446aa79634f99e8),
                },
                c1: Fp2 {
                    c0: Fp(0x1bb3e27460017785),
                    c1: Fp(0x2a043ae60a32e00b),
                },
                c2: Fp2 {
                    c0: Fp(0xd04fdbf60e75c4f),
                    c1: Fp(0x17b6ad9c501b4f3f),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xd83f1fcdcb18271),
                    c1: Fp(0x1317d762dd2f9ad7),
                },
                c1: Fp2 {
                    c0: Fp(0x1aba4a7167b7c8ed),
                    c1: Fp(0x2fa14821e0285b83),
                },
                c2: Fp2 {
                    c0: Fp(0x22e8380a4fa6cc60),
                    c1: Fp(0x229ad2dda454962b),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x138cf297ed68230f),
                    c1: Fp(0x1c7f0ba9ac7e4934),
                },
                c1: Fp2 {
                    c0: Fp(0x2c5582edac3c0f16),
                    c1: Fp(0xf0e8107881b6a38),
                },
                c2: Fp2 {
                    c0: Fp(0x9755169e28e6e5d),
                    c1: Fp(0xc1a82cca2b6eafe),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x3b406e5eeb3279),
                    c1: Fp(0x150091280aad7885),
                },
                c1: Fp2 {
                    c0: Fp(0xfa31fb2b23b3b96),
                    c1: Fp(0x24870dee3b8eee8c),
                },
                c2: Fp2 {
                    c0: Fp(0x220e2f63516ea26a),
                    c1: Fp(0x1ed289c2421692ca),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x9af8a3b28b97cd),
                    c1: Fp(0xa56d57961d72f7),
                },
                c1: Fp2 {
                    c0: Fp(0x23890112c1547343),
                    c1: Fp(0xbddfd15855b965b),
                },
                c2: Fp2 {
                    c0: Fp(0xe9dc6a17b1a7c2c),
                    c1: Fp(0x212fe840d3f61a7),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x2076f2ece91659f),
                    c1: Fp(0x24cb99d3acd38df4),
                },
                c1: Fp2 {
                    c0: Fp(0x14b49e0bf3d0ad20),
                    c1: Fp(0x37978fb25be2da7a),
                },
                c2: Fp2 {
                    c0: Fp(0x2ceb2d89f863904a),
                    c1: Fp(0xed8ff05ccd99da9),
                },
            },
            z: Fp6::one(),
        },
        ProjectivePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x31b620af57a5ebc8),
                    c1: Fp(0x3a7d3e886d2f9f59),
                },
                c1: Fp2 {
                    c0: Fp(0x340598d9017125),
                    c1: Fp(0x71afcaf8c6796cc),
                },
                c2: Fp2 {
                    c0: Fp(0xafbdcfde0edc5ac),
                    c1: Fp(0x263f8a4fc0cf0c98),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xae49a25c76c3535),
                    c1: Fp(0x352c8e59105f129e),
                },
                c1: Fp2 {
                    c0: Fp(0x1fb8ddda1167cc07),
                    c1: Fp(0xd76d44911e0763f),
                },
                c2: Fp2 {
                    c0: Fp(0x1d8ffb4fc31f2153),
                    c1: Fp(0x3672b55ffcbdbaf9),
                },
            },
            z: Fp6::one(),
        },
    ]),
]);
