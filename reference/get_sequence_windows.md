# extract sequence windows, from sequence at pos in protein

extract sequence windows, from sequence at pos in protein

## Usage

``` r
get_sequence_windows(
  unique_prot_pep_seq,
  flank_size = 7,
  sequence = "sequence",
  pos_in_protein = "pos_in_protein"
)
```

## Arguments

- unique_prot_pep_seq:

  data.frame with sequence and pos in protein

- flank_size:

  size of the window to extract

- sequence:

  name of column with sequence

- pos_in_protein:

  name of column with position in protein

## Examples

``` r
# Create sample data
sample_data <- data.frame(
  protein_id = c("P12345", "P12345", "Q67890"),
  sequence = c("MKFLVLLFNILCLFPVLAADNH", "MKFLVLLFNILCLFPVLAADNH", "AEQKLISEEDLLRKRREQLKHKLEQL"),
  pos_in_protein = c(5, 12, 8),
  peptide = c("FLV", "ILC", "EED")
)

# Extract sequence windows with default flank size (7)
result <- get_sequence_windows(sample_data)
expected <- c("XXXMKFLVLLFNILC", "VLLFNILCLFPVLAA", "AEQKLISEEDLLRKR")
stopifnot(all(result$sequence_window == expected))

# Extract sequence windows with custom flank size
result_small <- get_sequence_windows(sample_data, flank_size = 3)
stopifnot(all(result_small$sequence_window == c( "KFLVLLF",  "NILCLFP",  "LISEEDL")))

# Extract sequence windows with different column names
sample_data2 <- data.frame(
  prot_seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRIY",
  position = c(10, 25, 40)
)

result2 <- get_sequence_windows(sample_data2,
                               flank_size = 5,
                               sequence = "prot_seq",
                               pos_in_protein = "position")
print(result2$sequence_window)
#> [1] "YIAKQRQISFV" "SRQLEERLGLI" "PILSRVGDGTQ"
```
