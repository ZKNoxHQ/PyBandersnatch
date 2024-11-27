F = Field(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
a = F(26970114609044244057035295480476738649110987753426044830880824247613503541996)
d = F(5995764538993767865256199277202352314034766753214989701839360767638071068190)
r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
h = 4
E = CurveEdwards(a, d, r, h)
test_vectors = {}
test_vectors['p'] = E(F(0x3f7d7b15e4f4cf90518830d91ec62ff0610b6b5ff787f5320b7667e9dc79a6dc), F(
    0x70a5894a64445438d015ac32ba360f092cde44bab11fc2b7d4b5c0d216228cce), F(1))
test_vectors['q'] = E(F(0x4bdf65591a0db8728b543f87ffd3ab0ef8a6762cedc621fa1d7d9e51926d97e), F(
    0x6765f125acbdd10815bedbffe1147d36cf5c501f7ad467f4cd12aeec62a7ef1f), F(1))
test_vectors['p_double'] = E(F(0x3e0a2f142e41694880e6dbf17d4cc1e2592612a4992c272aec9b3f57f3d13834), F(
    0x3e8ff40606be05404b97ccdefb19956b75ab0b532da07edf5d6b4f6a404b06cf), F(1))
test_vectors['p_plus_q'] = E(F(0x299e56fa710434f85ba535e66de5520927af55434e6c01f882a5fe634a337871), F(
    0x2d04373af410d24b4bb950f8b3096928852705ab12c7e7b46611fac016328df9), F(1))
test_vectors['p_minus_q'] = E(F(0x7319225e4edaec81071418b47d1ba693967b6692894adef6240e784df7e4d654), F(
    0x31c47086ea9419d29477f86ad322c6dd0fa263febb3446cfb045aec50db059a3), F(1))
test_vectors['φ_p'] = E(F(0xe2585aca336c58ba1937a82601d7d2256868b44f3bddb764a1f841ad79da206), F(
    0x15ef05ebd97664593eb98170626b11599b3a5b0d0002b19a30dd389e78392579), F(1))
test_vectors['k_times_p'] = E(F(0x2ecefe5d96bb1eedcd7d44b963c22d47f5ccd90c2e9c05935c220123ace5f8a), F(
    0xf2693e9239ee3709661fbf6c908de99ce7a41f149cefebe5ef6fc2c292bb4c6), F(1))
test_vectors['k1_times_p'] = E(F(0x703e04d247e473ab4271e3091eb33a26bfa5f0ed30775c47b5984c3bfc0cc83e), F(
    0x1ea4e7674cc0fe71301e2d695d07ebcb64aa7f36a467c37edb74cd487a0b5681), F(1))
test_vectors['k2_times_p'] = E(F(0x2d530cbbae85b9d2dd500b69ec18da5768e89a80755289540a48fcb3ce9ba1c7), F(
    0x6831d061687693bae155eb19b5cc0a9694d431afcf8d87673796eb8066b35b6e), F(1))
test_vectors['k1_times_p_plus_k2_times_q'] = E(F(0x5611d5dd5ef0d41c33c069c2216a43cabd46374898be55891518d28368768d1c), F(
    0x28a725ede6d5ca8746af9dd277ebfef3ed1af24ee04036cef883043703df4b46), F(1))
test_vectors['k'] = 0x1a862619b8224e61eb24bb583c84ce04913064d37308623924c7a64fcdc9f191
test_vectors['k1'] = -0x2286ed83a0b1545d1b7788921e40bb14
test_vectors['k2'] = 0x6886451b4aa55294c626bb34d42e242
test_vectors['λ'] = 0x13b4f3dc4a39a493edf849562b38c72bcfc49db970a5056ed13d21408783df05
