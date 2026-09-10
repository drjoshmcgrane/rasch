# DIF with several person factors and repeated measures

``` r

library(rasch)
```

## Design the analysis around the person

Differential item functioning is a violation of the Rasch model’s
requirement of invariant comparison: item locations should not depend on
which persons respond (Rasch 1961). When several person factors are
relevant, they should enter one model. Repeated observations require a
further distinction: group is a between-person factor, whereas occasion
varies within person.

For a factor \\G\\, class interval \\C\\, and standardised residual
\\z\\, the single-factor model is

\\ z=\mu+G+C+G\mathbin{:}C+\varepsilon. \\

The \\G\\ term tests uniform DIF. The \\G\mathbin{:}C\\ term tests
non-uniform DIF. With several factors, `dif_anova` fits their terms
jointly and uses Type II sums of squares.

The following dataset has two observations per person. Item I03 has a
group shift, item I06 has an occasion shift, and item I05 shifts for
group B at the second occasion only, a group-by-occasion interaction.

``` r

set.seed(21)
N <- 320
difficulty <- seq(-1.5, 1.5, length.out = 8)
theta <- rnorm(N)
group <- rep(c("A", "B"), each = N / 2)

make_wave <- function(occasion_shift, interaction_shift) {
  shift <- matrix(0, N, 8)
  shift[group == "B", 3] <- 1.2
  shift[, 6] <- occasion_shift
  shift[group == "B", 5] <- interaction_shift
  matrix(rbinom(N * 8, 1,
                plogis(outer(theta, difficulty, "-") - shift)), N, 8)
}

X <- rbind(make_wave(0, 0), make_wave(1.0, 2.0))
colnames(X) <- sprintf("I%02d", 1:8)
dat <- data.frame(
  pid = rep(sprintf("P%03d", seq_len(N)), 2),
  X,
  group = rep(group, 2),
  occasion = rep(c("T1", "T2"), each = N)
)

# the same three persons at each occasion: the identifier repeats down the
# rows, which is what makes the design repeated measures
dat[c(1:3, N + 1:3), c("pid", "I01", "I02", "I03", "group", "occasion")]
#>      pid I01 I02 I03 group occasion
#> 1   P001   1   1   0     A       T1
#> 2   P002   1   1   1     A       T1
#> 3   P003   1   1   1     A       T1
#> 321 P001   0   1   0     A       T2
#> 322 P002   1   1   0     A       T2
#> 323 P003   1   1   1     A       T2
```

## Fit once and test both factors

The repeated person identifier is carried into the fit. `dif_anova`
detects occasion as within-person; it can also be declared explicitly.
Persons, not stacked rows, are the units of analysis. Uniform
between-person terms use Type II tests with HC3 covariance.
Class-interval interactions retain the residual-ANOVA reference.
Within-person tests use person-level contrasts, with a
Greenhouse–Geisser correction when a factor has more than two levels.
This mixed-design analysis extends the single-factor residual analysis
of variance described by Andrich and Marais (2019). Its F references are
large-sample approximations.

With incomplete panels, between-person tests fit the occasion and person
factors jointly, giving each person total weight one. They use
person-cluster CR3 covariance for both uniform and non-uniform terms,
including uncertainty in the occasion adjustment. This avoids confusing
unequal occasion coverage with differences in the composition of the
person groups. Within-person contrasts still require complete panels.

``` r

fit <- rasch(dat, id = "pid", factors = c("group", "occasion"),
             items = sprintf("I%02d", 1:8))
da <- dif_anova(fit, within = "occasion", effects = "factorial", sizes = TRUE)
da$summary
#>  item           term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#>   I01          group     0.786     0.376         1.000        0.003            
#>   I01       occasion     0.875     0.350         1.000        0.003            
#>   I01 group:occasion     3.153     0.077         1.000        0.010            
#>   I02          group     0.287     0.593         1.000        0.001            
#>   I02       occasion     4.885     0.028         1.000        0.016            
#>   I02 group:occasion     0.809     0.369         1.000        0.003            
#>   I03          group    15.018   < 0.001         0.006        0.048           *
#>   I03       occasion     1.340     0.248         1.000        0.004            
#>   I03 group:occasion     2.787     0.096         1.000        0.009            
#>   I04          group     1.292     0.257         1.000        0.004            
#>   I04       occasion     1.327     0.250         1.000        0.004            
#>   I04 group:occasion     1.821     0.178         1.000        0.006            
#>   I05          group    15.045   < 0.001         0.006        0.051           *
#>   I05       occasion    10.604     0.001         0.055        0.033            
#>   I05 group:occasion    27.711   < 0.001       < 0.001        0.083           *
#>   I06          group     6.249     0.013         0.544        0.020            
#>   I06       occasion    19.306   < 0.001       < 0.001        0.059           *
#>   I06 group:occasion     2.340     0.127         1.000        0.008            
#>   I07          group     4.837     0.029         1.000        0.015            
#>   I07       occasion     4.613     0.033         1.000        0.015            
#>   I07 group:occasion     2.761     0.098         1.000        0.009            
#>   I08          group     1.259     0.263         1.000        0.004            
#>   I08       occasion     0.491     0.484         1.000        0.002            
#>   I08 group:occasion     1.593     0.208         1.000        0.005            
#>  F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#>         0.212        0.932            1.000           0.003               
#>         1.165        0.326            1.000           0.015               
#>         1.051        0.381            1.000           0.014               
#>         3.400        0.010            0.416           0.042               
#>         0.814        0.517            1.000           0.010               
#>         1.631        0.166            1.000           0.021               
#>         2.355        0.054            1.000           0.030               
#>         2.148        0.075            1.000           0.027               
#>         0.421        0.793            1.000           0.005               
#>         0.513        0.726            1.000           0.007               
#>         1.739        0.141            1.000           0.022               
#>         1.083        0.365            1.000           0.014               
#>         2.483        0.044            1.000           0.031               
#>         0.812        0.518            1.000           0.010               
#>         2.916        0.022            0.886           0.037               
#>         1.294        0.272            1.000           0.017               
#>         1.333        0.257            1.000           0.017               
#>         1.000        0.408            1.000           0.013               
#>         0.558        0.693            1.000           0.007               
#>         0.645        0.631            1.000           0.008               
#>         1.597        0.175            1.000           0.020               
#>         1.514        0.198            1.000           0.019               
#>         0.386        0.818            1.000           0.005               
#>         0.886        0.472            1.000           0.011               
#>  superseded
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>           *
#>           *
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#> 
```

The multiplicity adjustment covers the complete family of
item-by-DIF-term tests. Uniform DIF is a factor effect that is stable
over the trait; a factor-by-class-interval effect is non-uniform DIF.
`effects = "factorial"` adds the person-factor interactions. A
significant higher-order term supersedes its component group terms
within the same item: item I05’s group effect is significant on its own,
but the `superseded` flag records that the interaction absorbs it, and
the follow-ups report the interaction rather than its components. Read
adjusted probabilities with effect sizes before changing an item.

[`dif_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
repeats this design under the fitted invariant model. Ordinary and
explanatory Rasch models and Multiple Ratings condition on each person’s
score. Extended Frames conditions within each item set, where the frame
unit is common. Comparative Judgement instead draws outcomes from its
fitted model while retaining judges and the comparison design. Every
replicate repeats the calibration and the full DIF family. The result is
a sensitivity analysis; the adjusted residual ANOVA remains primary. Its
familywise probabilities describe the fitted global invariant null, not
otherwise invariant items when another item has DIF.

``` r

db <- dif_bootstrap(fit, da, B = 999, workers = 4, seed = 2026)
db$summary[, c("item", "term", "p_uniform_boot_adj",
               "p_nonuniform_boot_adj")]
```

In a four-occasion validation study, familywise error under the balanced
Rasch null was 6.9% for the primary analysis and 4.9% for the bootstrap
over 350 datasets. For uniform occasion DIF the bootstrap was slightly
less powerful than the primary analysis. Adjusted non-uniform power was
weak for both procedures in this seven-item, 180-person design. The
bootstrap therefore remains a sensitivity analysis rather than the
default test.

The earlier panel-loss results used the superseded marginal adjustment.
For the current joint adjustment, a known-residual study with correlated
person factors and unequal panels gave uniform and non-uniform null
rejection rates of 5.7% and 5.3% over 3,000 replicates (Monte Carlo SEs
0.42 and 0.41 percentage points). A separate end-to-end Rasch check gave
3.0% over 200 datasets for a null group effect on an item carrying DIF
for another factor. These checks assess the primary tests, not the full
bootstrap family’s error rate under every incomplete-panel design.

## Quantify the departure

ANOVA identifies evidence against invariance; it does not state the size
of the departure in logits. `dif_size` resolves an item by a
between-person factor. `dif_contrasts` provides planned contrasts and
uses person-level differencing for within-person questions.

``` r

dif_size(fit, "I03", by = "group")
#> DIF size for I03 by group (resolved locations, logits)
#>  level location    se weak   n
#>      A   -0.518 0.129    0 160
#>      B    0.227 0.128    0 160
#>  level_a level_b difference    se      t  df       p   p_adj  lower  upper
#>        A       B     -0.745 0.203 -3.667 316 < 0.001 < 0.001 -1.145 -0.345
#>  significant practical p_beyond_A p_beyond_A_adj ets signed_area
#>            *   >= 0.50      0.058          0.058  B-            
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
dc <- dif_contrasts(fit, items = c("I03", "I06"), within = "occasion")
dc$table
#>  item                         contrast within estimate se statistic      df
#>   I03                     group: B - A           0.737        4.044 316.916
#>   I03                occasion: T2 - T1      *   -0.258       -0.876 313.935
#>   I03 group(B - A) x occasion(T2 - T1)      *   -0.258       -1.397 313.935
#>   I06                     group: B - A          -0.434       -2.138 308.822
#>   I06                occasion: T2 - T1      *    0.906        4.438 298.364
#>   I06 group(B - A) x occasion(T2 - T1)      *    0.134        1.220 298.364
#>        p   p_adj lower upper significant practical
#>  < 0.001 < 0.001                       *         *
#>    0.382   0.490                                  
#>    0.163   0.490                                  
#>    0.033   0.133                                  
#>  < 0.001 < 0.001                       *         *
#>    0.223   0.490
da$posthoc
#>  item           term item.1        contrast within estimate se statistic
#>   I03          group    I03           B - A           0.737        4.044
#>   I05 group:occasion    I05 B - A x T2 - T1      *    2.140        5.137
#>   I06       occasion    I06         T2 - T1      *    0.906        4.438
#>       df       p   p_adj lower upper significant practical
#>  316.916 < 0.001 < 0.001                       *         *
#>  317.916 < 0.001 < 0.001                       *         *
#>  298.364 < 0.001 < 0.001                       *         *
```

[`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
is the general follow-up for a significant term. A main effect with more
than two levels is reported as pairwise marginal differences over the
other fitted factors. An interaction is reported as a
difference-in-differences, or its higher-order counterpart. These
comparisons use the joint covariance of the resolved item locations and
Holm adjustment over the stated family. The result is on the logit scale
and respects the factor structure used in the DIF analysis. Here item
I03 is reported as a pairwise group difference, item I06 as a
person-level occasion contrast, and item I05 as the
difference-in-differences of its interaction; the superseded I05 group
term receives no follow-up of its own.

For repeated-person contrasts, significance comes from person-level
residual contrast scores with the same design-cell weights as the
resolved estimate. Other fitted factors are averaged equally over their
complete cells, including when their sample sizes differ, and the
independent between-person cells use a Welch–Satterthwaite reference.
The resolved logit difference remains the magnitude. Its standard error
uses the person-clustered calibration covariance, whereas
[`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
tests a nominated within-person comparison directly from person-level
contrasts. The fitted person identifier is used unless `id` is supplied
explicitly.

For a many-facet fit, follow-ups may name the underlying item; its
virtual facet cells are pooled with common weights so facet severity
cancels from the group contrast. Extended-frame fits support the
residual ANOVA for factors outside the frame definition, but not an
ordinary resolved-item magnitude: that refit would discard the fitted
frame units.

A statistical flag should be considered with the logit magnitude,
targeting, item content, and the intended use of the scale. Resolving an
item changes the measurement model and should follow a substantive
account of why the item is not invariant.

## References

Andrich, D., and Marais, I. (2019). *A Course in Rasch Measurement
Theory: Measuring in the Educational, Social and Health Sciences*.
Springer.

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In *Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 4, pp. 321–333).
Berkeley: University of California Press.
