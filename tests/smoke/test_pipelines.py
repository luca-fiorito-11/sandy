def test_xs_pipeline():
    import sandy
    tape = sandy.get_endf6_file("jeff_33", "xs", 10010).to_file("H1.jeff33")
    cli = "H1.jeff33 --samples 2  --acer --temperatures 300".split()
    sandy.sampling.run(cli)


def test_fy_pipeline():
    import sandy
    tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350).to_file("U5.jeff33")
    cli = "U5.jeff33 --samples 2  --fycov  --mt 102".split()
    sandy.sampling.run(cli)


def test_rdd_pipeline():
    import sandy
    tape = sandy.get_endf6_file("jeff_33", "decay", 551340).to_file("Cs134.jeff33")
    cli = "Cs134.jeff33 --samples 2".split()
    sandy.sampling.run(cli)