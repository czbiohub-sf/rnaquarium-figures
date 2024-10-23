from scripts.util.fuzzy import classify_stage

#########
# TESTS
#########

def test_hfp_a():
    text = "Example with hfp 2.5"
    n, t = classify_stage(text)
    assert (n == 2.5) and (t == "hfp")

def test_dfp_a():
    text = "Example with dfp 2.5"
    n, t = classify_stage(text)
    assert (n == 2.5) and (t == "dfp")
    
def test_dfp_b_c_1():
    text = "Example with 0.75mm length but also with 2.5 dfp"
    n, t = classify_stage(text)
    assert (n == 2.5) and (t == "dfp")

def test_dfp_b_c_2():
    text = "Example with 0.75mm length but also with 1 dfp"
    n, t = classify_stage(text)
    assert (n == 1) and (t == "dfp")

def test_dfp_b_c_3():
    text = "Example with 0.75mm length but also with 1.5 dfp"
    n, t = classify_stage(text)
    assert (n == 1.5) and (t == "dfp")

def test_hpf_b_c_1():
    text = "Example with 0.75mm length but also with 1.5 hpf"
    n, t = classify_stage(text)
    assert (n == 1.5) and (t == "hpf")

def test_hour_b_c_1():
    text = "Example with 0.75mm length but also with 1.5 hour"
    n, t = classify_stage(text)
    assert (n == 1.5) and (t == "hour")

def test_somite_b_c_1():
    text = "Example with 0.75mm length but also with 15 somites"
    n, t = classify_stage(text)
    assert (n == 15) and (t == "somites")

def test_somite_b_c_2():
    text = "Example with 0.75mm length but 4 somite stage"
    n, t = classify_stage(text)
    assert (n == 4) and (t == "somite")

def test_hours_fert_1():
    text = "collected at 48 hours post-fertlilization"
    n, t = classify_stage(text)
    assert (n == 48) and (t == "hours post-fertlilization")

def test_days_fert_2():
    text = "mutant zebrafish larvae at 3 days post fertilisation were pooled into three biological replicates each of 50 larvae."
    n, t = classify_stage(text)
    assert (n == 3) and (t == "days post fertilisation")

def test_days_fert_3():
    text = "RNA-Sequencing analysis of 3 days post fertlisation whole wild-type and prp1-/-;prp2-/- zebrafish larvae"
    n, t = classify_stage(text)
    assert (n == 3) and (t == "days post fertlisation")

def test_hours_fert_4():
    text = "digestion of whole embryos at 48 hours post-fertlization."
    n, t = classify_stage(text)
    assert (n == 48) and (t == "hours post-fertlization")