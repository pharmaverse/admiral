
derive_var_eosdy <- function(dataset, start_date = TRTSDT, end_date = EOSDT) {
  derive_duration(
    dataset,
    new_var = EOSDY,
    start_date = !!enquo(start_date),
    end_date = !!enquo(end_date)
  )
}
