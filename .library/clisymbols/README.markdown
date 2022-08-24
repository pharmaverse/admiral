



# clisymbols

[![Linux Build Status](https://travis-ci.org/gaborcsardi/clisymbols.svg?branch=master)](https://travis-ci.org/gaborcsardi/clisymbols)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/github/gaborcsardi/clisymbols?svg=true)](https://ci.appveyor.com/project/gaborcsardi/clisymbols)
[![](http://www.r-pkg.org/badges/version/clisymbols)](http://www.r-pkg.org/pkg/clisymbols)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/clisymbols)](http://www.r-pkg.org/pkg/clisymbols)


> Unicode symbols with Windows fallbacks

![](/inst/screenshot.png)

Inspired by (and mostly copied from) the
[figures](https://github.com/sindresorhus/figures) JavaScript project.

## Install


```r
install.packages("devtools")
devtools::install_github("gaborcsardi/clisymbols")
```

## Usage


```r
library(clisymbols)
cat(symbol$tick, "All good\n")
```

```
#> ✔ All good
```

```r
cat(symbol$cross, "Problem\n")
```

```
#> ✖ Problem
```

Here is a list of all symbols, with their names:


```r
for (i in seq_along(symbol)) {
  cat(symbol[[i]], "\t", names(symbol)[i], "\n", sep = "")
}
```

```
#> ✔	tick
#> ✖	cross
#> ★	star
#> ▇	square
#> ◻	square_small
#> ◼	square_small_filled
#> ◯	circle
#> ◉	circle_filled
#> ◌	circle_dotted
#> ◎	circle_double
#> ⓞ	circle_circle
#> ⓧ	circle_cross
#> Ⓘ	circle_pipe
#> ?⃝	circle_question_mark
#> ●	bullet
#> ․	dot
#> ─	line
#> ═	double_line
#> …	ellipsis
#> ❯	pointer
#> ℹ	info
#> ⚠	warning
#> ☰	menu
#> ☺	smiley
#> ෴	mustache
#> ♥	heart
#> ↑	arrow_up
#> ↓	arrow_down
#> ←	arrow_left
#> →	arrow_right
#> ◉	radio_on
#> ◯	radio_off
#> ☒	checkbox_on
#> ☐	checkbox_off
#> ⓧ	checkbox_circle_on
#> Ⓘ	checkbox_circle_off
#> ❓	fancy_question_mark
#> ≠	neq
#> ≥	geq
#> ≤	leq
#> ▔	upper_block_1
#> ▀	upper_block_4
#> ▁	lower_block_1
#> ▂	lower_block_2
#> ▃	lower_block_3
#> ▄	lower_block_4
#> ▅	lower_block_5
#> ▆	lower_block_6
#> ▇	lower_block_7
#> █	lower_block_8
#> █	full_block
```

### Fallback symbols

Some terminals do not support (all) Unicode characters, and on these reasonable
ASCII substitutes are used:


```
#> √   tick
#> x   cross
#> *   star
#> █   square
#> [ ] square_small
#> [█] square_small_filled
#> ( ) circle
#> (*) circle_filled
#> ( ) circle_dotted
#> (o) circle_double
#> (o) circle_circle
#> (x) circle_cross
#> (|) circle_pipe
#> (?) circle_question_mark
#> *   bullet
#> .   dot
#> ─   line
#> =   double_line
#> ... ellipsis
#> >   pointer
#> i   info
#> ‼   warning
#> ≡   menu
#> ☺   smiley
#> ┌─┐ mustache
#> ♥   heart
#> ^   arrow_up
#> v   arrow_down
#> <   arrow_left
#> >   arrow_right
#> (*) radio_on
#> ( ) radio_off
#> [x] checkbox_on
#> [ ] checkbox_off
#> (x) checkbox_circle_on
#> ( ) checkbox_circle_off
#> (?) fancy_question_mark
#> !=  neq
#> >=  geq
#> <=  leq
#> ^   upper_block_1
#> ^   upper_block_4
#> .   lower_block_1
#> _   lower_block_2
#> _   lower_block_3
#> =   lower_block_4
#> =   lower_block_5
#> *   lower_block_6
#> █   lower_block_7
#> █   lower_block_8
#> █   full_block
```

# License

MIT © [Gabor Csardi](http://gaborcsardi.org) and [Sindre Sorhus](http://sindresorhus.com)
