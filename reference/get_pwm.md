# Create PWM from amino acid sequences

Internal helper function to compute a position weight matrix (PWM) from
a vector of amino acid sequences.

## Usage

``` r
get_pwm(sequences)
```

## Arguments

- sequences:

  Character vector of amino acid sequences (must be same length)

## Value

A matrix with amino acids as rows and positions as columns
