import phydrus as ps


def test_read_profile():
    ps.read_profile(path="tests/test_data/PROFILE.OUT")
    return


def test_read_nod_inf():
    ps.read_nod_inf(path="tests/test_data/NOD_INF.OUT")
    return


def test_read_run_inf():
    ps.read_run_inf(path="tests/test_data/RUN_INF.OUT")
    return


def test_read_obs_node():
    ps.read_obs_node(path="tests/test_data/OBS_NODE.OUT", nodes=[5, 7, 9])
    return


def test_read_i_check():
    ps.read_i_check(path="tests/test_data/I_CHECK.OUT")
    return


def test_read_tlevel():
    ps.read_tlevel(path="tests/test_data/T_LEVEL.OUT")
    return


def test_read_alevel():
    ps.read_alevel(path="tests/test_data/A_LEVEL.OUT")
    return


def test_read_balance():
    ps.read_balance(path="tests/test_data/BALANCE.OUT")
    return


def test_read_solute():
    ps.read_solute(path="tests/test_data/SOLUTE1.OUT")
    return
