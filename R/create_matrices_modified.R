#' Expand compound constraints
#'
#' Takes a compound BAIN constraint, with multiple (in)equality operators, and
#' expands it into simple cosource('E:/SSDANOVA/create_matrices_modified.R')nstraints.
#' @param hyp Character. A BAIN (in)equality constraint
#' @return A character vector with one element for each simple constraint
#' @examples
#' expand_compound_constraints("a=b=c")
#' expand_compound_constraints("openness>neuroticism>extraversion==c")
#' @keywords internal
expand_compound_constraints <- function(hyp){
  equality_operators <- gregexpr("[=<>]", hyp)[[1]]
  if(length(equality_operators) > 1){
    string_positions <- c(0, equality_operators, nchar(hyp)+1)
    return(sapply(1:(length(string_positions)-2), function(pos){
      substring(hyp, (string_positions[pos]+1), (string_positions[pos+2]-1))
    }))
  } else {
    return(hyp)
  }
}

#' Expand parentheses
#'
#' Takes a BAIN constraint with parentheses containing a "vector" of parameter
#' labels, and recursively expands the parenthesized "vectors".
#' @param hyp Character. A BAIN (in)equality constraint
#' @return A character vector with one element for each simple constraint
#' @examples
#' expand_parentheses("(a,b)>c")
#' expand_parentheses("(openness, conscientiousness)>(neuroticism,extraversion)")
#' @keywords internal
expand_parentheses <- function(hyp){
  parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
  if(!parenth_locations[1] == -1){
    if(length(parenth_locations) %% 2 > 0) stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
    expanded_contents <- strsplit(substring(hyp, (parenth_locations[1]+1), (parenth_locations[2]-1)), ",")[[1]]
    if(length(parenth_locations) == 2){
      return(paste0(substring(hyp, 1, (parenth_locations[1]-1)), expanded_contents, substring(hyp, (parenth_locations[2]+1), nchar(hyp))))
    } else {
      return(apply(
        expand.grid(expanded_contents, expand_parentheses(substring(hyp, (parenth_locations[2]+1), nchar(hyp)))),
        1, paste, collapse = ""))
    }
  } else {
    return(hyp)
  }
}

#' Flip inequality
#'
#' Takes a BAIN constraint and flips elements on both sides of any "<" operators
#' so that only the ">" operator is used to define inequality constraints.
#' @param hyp Character. A BAIN (in)equality constraint
#' @return Character
#' @examples
#' flip_inequality("b<c")
#' @keywords internal
flip_inequality <- function(hyp){
  if(grepl("<", hyp)){
    loc <- gregexpr("<", hyp)[[1]][1]
    return(paste0(substring(hyp, (loc+1)), ">", substring(hyp, 1, (loc-1))))
  } else {
    return(hyp)
  }
}

#' Constraint to equation
#'
#' Formats a BAIN constraint as an equation that can be evaluated. Adds scalars
#' to all parameters in the constraint, and adds a XXXconstant parameter to any
#' constants. Also adds "*" and "+" operators where necessary for evaluation.
#' @param hyp Character. A BAIN (in)equality constraint
#' @return Character
#' @examples
#' constraint_to_equation(hyp = c("-strength>wisdom", "intelligence>wisdom",
#'                     "strength>dexterity", "intelligence>dexterity"))
#' constraint_to_equation("1*b")
#' constraint_to_equation(hyp  ="b")
#' constraint_to_equation(hyp = "5*b-4c")
#' constraint_to_equation("5+c=0")
#' constraint_to_equation("5+c=d")
#' @keywords internal
constraint_to_equation <- function(hyp){
  # When the string starts with a word, OR when a word is not preceded by
  # a number or *-sign, replace the "word" with "1*word"
  # hyp <- gsub("(^|(?<![\\*\\d]))([a-zA-Z][a-zA-Z0-9_]{0,})", "1*\\2", hyp, perl = TRUE)
  # Gu: add if not being preceded by any of .a-zA-Z0-9 , then ignore. otherwise x.A becomes x.1*A
  hyp <- gsub("(^|(?<![\\*\\d\\.a-zA-Z0-9_]))([a-zA-Z][a-zA-Z0-9_]{0,})", "1*\\2", hyp, perl = TRUE)

  # If a number starts with a period, add a leading zero
  # hyp <- gsub("(^|(?<!\\d))(\\.[0-9]+)", "0\\2", hyp, perl = TRUE)
  # Gu: add if being preceded by any of a-zA-Z0-9_, then ignore. Otherwise x.1 becomes x0.1
  hyp <- gsub("(^|(?<![a-zA-Z0-9_]))(\\.[0-9]+)", "0\\2", hyp, perl = TRUE)

  # When a number is directly followed by a word, insert a * sign
  # hyp <- gsub("(\\d)(?=[a-zA-Z][a-zA-Z0-9_]{0,})", "\\1*", hyp, perl = TRUE)
  # Gu: but not precededed by any of a-zA-Z_. Otherwise, x1x becomes x1*x
  hyp <- gsub("(?<![a-zA-Z_\\.])(\\d)(?=[a-zA-Z][a-zA-Z0-9_]{0,})", "\\1*", hyp, perl = TRUE) #####

  # When the string starts with a floating point number, OR
  # when a floating point number is not preceded by + or -, add a +
  # hyp <- gsub("(^|(?<![+-]))([0-9]{1,}\\.[0-9]{0,})", "+\\2", hyp, perl = TRUE)
  # Gu: add a condition: followed by a *, because otherwise x0.1x becomes x+0.1x.
  hyp <- gsub("(^|(?<![+-]))([0-9]{1,}\\.[0-9]{0,}(?=\\*))", "+\\2", hyp, perl = TRUE)

  # When the string starts with an integer number, OR
  # when an integer number is not preceded by + or -, add a +
  # hyp <- gsub("(^|(?<![\\.0-9+-]))([0-9]{1,})(?!\\.)", "+\\2", hyp, perl = TRUE)
  # Gu: add a condition: followed by a *, because otherwise x1x becomes x+1x.
  hyp <- gsub("(^|(?<![\\.0-9+-]))([0-9]{1,})(?!\\.)(?=\\*)", "+\\2", hyp, perl = TRUE)

  # Gu: To complement the above two, for a single number (constant), add a +
  hyp <- gsub("(?<=[=<>]|^)(\\d|(\\d+\\.\\d+))(?=[=+-<>]|$)", "+\\1", hyp, perl = TRUE)

  # When a number is followed by =, >, <, +, or -, or is the last character of
  # the string, add "*XXXconstant" to the number
  # hyp <- gsub("(\\d)((?=[=<>+-])|$)", "\\1*XXXconstant", hyp, perl = TRUE)
  # Gu: but not preceded by word, number, period and  underline
  gsub("(?<![a-zA-Z0-9_\\.])(\\d|(\\d+\\.\\d+))((?=[=<>+-])|$)", "\\1*XXXconstant", hyp, perl = TRUE)
}

#' Order terms
#'
#' Moves all terms on the right hand side of a BAIN constraint, which has been
#' formatted as an equation, to the left hand side.
#' @param hyp Character. A BAIN (in)equality constraint formatted as equation.
#' @return Character
#' @examples
#' order_terms("+2+1*a-1=+1*b-2")
#' order_terms(hyp = "+1*a=-2")
#' order_terms(hyp = "+1*a=-2+3")
#' order_terms(hyp = "+1*a=+1*b")
#' order_terms(hyp = "-2+1*a=+1*b")
#' @keywords internal
order_terms <- function(hyp){
  eq_location <- gregexpr("[=<>]", hyp)[[1]]
  rhs <- substring(hyp, eq_location+1, nchar(hyp))
  rhs <- gsub("(?=[+-])", "-", rhs, perl = TRUE)
  paste0(rhs, substring(hyp, 1, eq_location))
}

#' Constraint to row
#'
#' Evaluate a BAIN constraint, which has been formatted as an equation, with all
#' terms moved to the left hand side, and return a single row of a BAIN
#' (in)equality constraint matrix.
#' @param hyp Character. A BAIN (in)equality constraint formatted as equation,
#' with all terms on the left hand side.
#' @return Numeric vector.
#' @examples
#' constraint_to_row(varnames = c("a", "b", "c", "d"), hyp = "+1*a-1*b+0=")
#' constraint_to_row(c("a", "b", "c", "d"), "+1*b-1*c+0=")
#' constraint_to_row(c("a", "b", "c", "d"), "-+0+1*c-1*d=")
#' @keywords internal
constraint_to_row <- function(varnames, hyp){
  e <- new.env()
  objects <- c(varnames, "XXXconstant")
  invisible(sapply(objects, assign, value = 0, envir = e))
  constraint_expression <- parse(text = substring(hyp, 1, nchar(hyp)-1))
  e[["XXXconstant"]] <- -1
  equal_to <- eval(constraint_expression, envir = e)
  e[["XXXconstant"]] <- 0

  c(sapply(varnames, function(x){
    e[[x]] <- 1
    constant <- eval(constraint_expression, envir = e)
    e[[x]] <- 0
    constant
  }), equal_to)
}

#' Create BAIN (in)equality constraint matrices
#'
#' Parses a character string describing a set of BAIN informative hypotheses,
#' and returns BAIN (in)equality constraint matrices. See Details for more
#' information.
#' @details Informative hypotheses specified as a character string should adhere
#' to the following simple syntax:
#' \itemize{
#'   \item Competing hypotheses are separated by ";".
#'         Thus, "a=b;a>b" means that H1: a=b, and H2: a>b.
#'   \item Each individual hypothesis consists of a (series of) (in)equality
#'         constraint(s). Every single (in)equality constraint is of the form
#'         "R1*mu1 + R2*mu2+... = r", where capital Rs refer to numeric scaling
#'         constants, must refer to the names of parameters in the model, and
#'         the lower case r refers to a constant. Standard mathematical
#'         simplification rules apply; thus, "R1*mu1 = R2*mu2" is equivalent to
#'         "R1*mu1 - R2*mu2 = 0".
#'   \item Multiple unrelated constraints within one hypothesis can be chained
#'         by "&". Thus, "a=b&c=d" means that H1: a=b AND c=d.
#'   \item Multiple related constraints within one hypothesis can be chained by
#'         repeating the (in)equality operators "=", "<", or ">". Thus, "a<b<c"
#'         means that H1: a < b AND b < c.
#'   \item Parameters can be grouped by placing them in a parenthesized, comma
#'         separated list. Thus, "(a,b)>c" means that H1: a > c AND b > c.
#'         Similarly, "(a,b)>(c,d)" means that H1: a > c AND b > c AND b > c AND
#'         b > d.
#' }
#' @param object Object of class \code{lm}, from which the model parameters are
#' extracted.
#' @param hyp Character string, containing a BAIN hypothesis (see Details).
#' @return A pair of named matrices for every hypothesis specified in the
#' \code{hyp} argument; one matrix named ERr, specifying equality constraints,
#' and one matrix named IRr, specifying inequality constraints.
#' #@import stats
#' #@import ranger
#' #@import metafor
#' @export
#' @examples
#' varnames <- c("strength","intelligence","wisdom","dexterity","constitution","charisma")
#' hyp <- "(strength, intelligence) > (wisdom, dexterity)"
#' create_matrices(varnames, hyp)
#' hyp2 <- hyp <- "strength > charisma > intelligence = .5"
#' create_matrices(varnames, hyp2)
#' hyp3 <- "wisdom > 0& (strength, intelligence) > dexterity = constitution"
#' create_matrices(varnames, hyp3)
#'
#' varnames <- c("a","b","c","d","e","f")
#' hyp <- "e=f<a=b=c"
#' create_matrices(varnames, hyp)
#'
#' hyp <- ".5*f>e>0&c>0&d=0&a=b=0"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "f>e<0&c>0&d=0&a=b=0"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "f<e>0&c>0&d=0&a=b=0"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "a>b<c>d"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "d<c>b<a"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "a>2"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "a>2& b+3=0& b>c=d&e=f=4"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "(a,b)>c&d=0"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "a>-1&a<1"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "a>b&b>c&c>a"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "a+b>2"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "(a+b)>2"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "a&b>c&d"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "(.5*a,b)>c"
#' create_matrices(varnames, hyp)
#'
#' hyp <- "1/2*a>c"
#' create_matrices(varnames, hyp)
#'
#' varnames <- c("a", "b", "c")
#' hyp <- "a < -2"
#' create_matrices(varnames, hyp)
#' hyp <- "-2 > a"
#' create_matrices(varnames, hyp)
#'
#' # HET - TEKEN VOOR DE C CONSTRAINT WORD IN DE UITVOER NIET NAAR EEN + OMGEZET
#' varnames <- c("a","b","c","d","e","f")
#' hyp <- "a>2; b=0; c< -0.5; d>e=f"
#' create_matrices(varnames, hyp)
#'
#' # -.5 WORD NIET BEGREPEN, -0.5 WEL
#' hyp1 <- "a>2; b=0; c< -.5; d>e=f"
#' create_matrices(varnames, hyp1)
#'
#' varnames <- c("a","b","c","d","e","f")
#' hyp <- "(a+a+a,c+c+c)>b+b"
#' create_matrices(varnames, hyp)
#' varnames <- c("a","b","c","d","e","f")
#' hyp1 <- "a+a+b>b+b+a"
#' create_matrices(varnames, hyp1)
create_matrices <- function(varnames, hyp){
  if(is.null(varnames)) stop("Please input proper linear model object")
  hyp <- gsub("\\s", "", hyp)
  if(grepl("[><=]{2,}", hyp)) stop("Do not use combined comparison signs e.g., '>=' or '=='")

  hyp_list <- strsplit(hyp, ";")[[1]]
  hyp_list <- lapply(hyp_list, function(x){ strsplit(x, "&")[[1]]})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_compound_constraints))})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_parentheses))})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, flip_inequality)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, constraint_to_equation)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, order_terms)})

  hyp_list <- unlist(lapply(hyp_list, function(x){
    ERr <- x[grep("=", x)]
    IRr <- x[grep("[<>]", x)]
    if(length(ERr) == 0){
      ERr <- NULL
    } else {
      ERr <- t(sapply(ERr, constraint_to_row, varnames = varnames))
      colnames(ERr)[ncol(ERr)]<- "="
    }
    if(length(IRr) == 0){
      IRr <- NULL
    } else {
      IRr <- t(sapply(IRr, constraint_to_row, varnames = varnames))
      colnames(IRr)[ncol(IRr)]<- ">"
    }
    list(ERr = ERr, IRr = IRr)
  }), recursive = FALSE)

  names(hyp_list) <- paste0(names(hyp_list), rep(1:(length(hyp_list)/2), each = 2))
  hyp_list
}
