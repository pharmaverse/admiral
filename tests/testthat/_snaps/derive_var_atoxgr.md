# derive_var_atoxgr Test 129: CTCAEv6  Blood bilirubin increased

    Code
      expect_dfs_equal(base = expected_bili_ctcv6, compare = actual_bili_ctcv6, keys = c(
        "ATOXDSCH", "AVAL", "ANRHI", "BNRIND", "AVALU"))

