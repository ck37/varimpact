# Function that counts # of unique values.
# Do not export.
length_unique = function(x) {
  # Skip NAs.
  length(setdiff(unique(x), NA))
}
