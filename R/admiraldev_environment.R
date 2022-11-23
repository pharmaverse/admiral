admiraldev_environment <- new.env(parent = emptyenv())

# assertions.R
admiraldev_environment$many_to_one <- NULL
admiraldev_environment$one_to_many <- NULL

# datasets.R
# get_dataset() is used to retrieve many_to_one and one_to_many
