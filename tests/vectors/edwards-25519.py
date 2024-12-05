# File generated using `sage sage/edwards.sage > tests/vectors/edwards.py`.
from src.field import Field
from src.curve.edwards25519 import Edwards25519
F = Field(0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed)
a = F(57896044618658097711785492504343953926634992332820282019728792003956564819944)
d = F(11790395817372903580333940030740964167805592400755249022757551653560006957928)
r = 7237005577332262213973186563042994240857116359379907606001950938285454250989
h = 8
E = Edwards25519(a, d, r, h)
test_vectors = {}
test_vectors['p'] = E(F(0x1cafec96dcbc2a677909d15e954584c9b679626a5ea497b285ef3e16e6a81cb9), F(0x7dbbb699c789815f8e9072ff1d83b1943fb3276eed60dcbfe860af964dac8648), F(1))
test_vectors['q'] = E(F(0x6291b7e19fd605c03dabab272180f94f6f6c95c823dd93bc4b14b9f8e93c684b), F(0x47bb7abaf9b01f1dc156d5454800ceaec580dbf7826127b39d8a5608e4e217d2), F(1))
test_vectors['p_double'] = E(F(0x12fa2884427f8791cf2207cb86bc3cc9637b412d4baefa0ec8a5422216ecec3a), F(0x25a7e0fe3b4ecea3587d11ecffff1be0bda2c3c1f9a2eb497acc7ec15f4f682f), F(1))
test_vectors['p_plus_q'] = E(F(0x2ff0ef8ed84b0f3a0282cc87fde82783d744c0d9621997b21f4ef691db2fa06d), F(0x201a64d1a59981b16ecfddc16921c0fd0258abc6d46ce2ff21dd2c6aa6f4bcae), F(1))
test_vectors['p_minus_q'] = E(F(0x44d627e115da9805ddee7d54e4aae5db69b4137e4f624b7f7b87cacec54b9291), F(0x213ead4de38cdaf53b1355fddb781dd3188b972636b3af0e48c7b68822b46ff1), F(1))
test_vectors['k_times_p'] = E(F(0xe29761e1f6f46e0a673344963de99125c2fa216e479e2d7f995a35a686561f7), F(0x656bd7c0f0f89e2ca28e85b187b13b4a65a242d7b6d579b9c7b25f6a37a64c52), F(1))
test_vectors['k1_times_p'] = E(F(0xa26ebbf28bd9902cb4a3eb7d11294ec581a1de3de1c10bb458dbea203a73528), F(0x3fbf2886e78dd069e55cc88d0cbd29bd4f029fd06d742809299265447b7dc6ae), F(1))
test_vectors['k2_times_p'] = E(F(0x443793f05874d25dd1d1999f0ed83f170ac244a146dc84b682f014e20c3fc29c), F(0x498fb795c5b7c0e942b2ff0e1c558fd78b6724881cc66a976661cbd9f5289063), F(1))
test_vectors['k1_times_p_plus_k2_times_q'] = E(F(0x5da846ad6c7c0c4e3683aeaf6940814e472efc2bb78b1762f245aa8e636ce052), F(0x2e2019dbce08b0491311162d601eb53adbf5eb9433c2260c427636736ad0605e), F(1))
test_vectors['k'] = 0x1a862619b8224e61eb24bb583c84ce04913064d37308623924c7a64fcdc9f191
test_vectors['k1'] = -0x2286ed83a0b1545d1b7788921e40bb14
test_vectors['k2'] = 0x6886451b4aa55294c626bb34d42e242
test_vectors['p_order_2_1'] = E(F(0x0), F(0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffec), F(1))
test_vectors['p_order_2_2'] = E(F(0x0), F(0x1), F(0x0))
test_vectors['p_order_2_3'] = E(F(0x1), F(0x0), F(0x0))
