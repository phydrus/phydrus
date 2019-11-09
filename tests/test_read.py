import pydrus as ps


def test_read_profile():
    ps.read_profile(path="test_data/PROFILE.OUT")
    return


def test_read_nod_inf():
    ps.read_nod_inf(path="test_data/NOD_INF.OUT")
    return


def test_read_run_inf():
    ps.read_run_inf(path="test_data/RUN_INF.OUT")
    return


def test_read_obs_node():
    ps.read_obs_node(path="test_data/OBS_NODE.OUT", nodes=[5, 7, 9])
    return


def test_read_i_check():
    ps.read_i_check(path="test_data/I_CHECK.OUT")
    return


def test_read_tlevel():
    ps.read_tlevel(path="test_data/T_LEVEL.OUT")
    return


def test_read_alevel():
    ps.read_alevel(path="test_data/A_LEVEL.OUT")
    return


def test_read_balance():
    ps.read_balance(path="test_data/BALANCE.OUT")
    return


def test_read_solute():
    ps.read_solute(path="test_data/SOLUTE1.OUT")
    return
