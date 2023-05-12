# Functions

#' @title Taylor Remainder
#' @description
#' The Taylor Remainder Function is integrated with the make_taylor_series
#' function to provide an error estimate when calculating the Taylor series
#' for a given function.
#'
#' @param df_np1 The n+1th derivative of the last derivative present in the
#' Taylor polynomial
#' @param x the value of x to plug into the absolute value of the
#' n+1th derivative
#' @param n The integer number of derivatives taken in total
#' @param a The center being tested
#' @returns Returns the absolute value of n+1th derivative multiplied by
#' x - a raised to a power of n+1 all divided by the factorial of n+1
#' @export
R_nx <-
  \(df_np1, x, n, a) {
    (abs(eval(df_np1, list(x = x))) / factorial(n)) * (abs(x - a) ^ (n))
  }

#' @title Taylor Series Converter
#' @description
#' The Taylor Series Converter Function performs a single variable Taylor
#' Expansion of the given function. It returns the call in full, the error
#' approximation, and the Taylor Remainder graph.
#'
#' @param f The function you want to convert into a Taylor Series
#' @param num_derivs How derivatives you want to take of f
#' @param center The center to test. In Stewart Calculus, this is notated "a"
#' @param intr_start The lower bound to try for the remainder interval
#' @param intr_end The upper bound to try for the remainder interval
#' @param var_name Differentiating f w.r.t var_name. Can only be one variable name
#' at the moment. The default variable name is "x".
#' @param open_closed Whether or not the interval defined by sequence
#' from intr_start to intr_end
#' @export
make_taylor_series <-
  \(
    f,
    num_derivs,
    center,
    intr_start,
    intr_end,
    var_name = "x",
    open_closed = "closed"
  ) {
    # Redefine terms
    a <- center
    c <- intr_start
    d <- intr_end

    # Parse Body of supplied function
    exprn <-
      body(f) %>%
      as.character() %>%
      .[[2]] %>%
      parse(text = .)

    # Run the symbolic differentiation
    d_list <- list()
    for (i in seq_along(1:(num_derivs + 1))) {
      if (i == 1) {
        d_list[[i]] <- exprn %>% .[[1]]
      } else if (i == 2) {
        d_list[[i]] <- exprn %>% stats::D(var_name)
      } else {
        d_list[[i]] <- d_list[[i - 1]] %>% stats::D(var_name)
      }
    }
    test <- f(a)
    center_check <- is.na(test) | is.infinite(test) | is.nan(test)
    if (center_check) {
      a <- a + 1e-6
    }
    # Calculate the terms for the Taylor series
    call_tbl <-
      tibble::tibble(n = 0:(length(d_list) - 1)) %>%
      dplyr::mutate(a = a, call = d_list) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(term_1 = eval(call, list(x = a)) / factorial(n)) %>%
      dplyr::mutate(
        term_2 =
          glue::glue("(", "(x-a)^", "{n})") %>%
          stringr::str_replace("\\ba\\b", as.character(a)),
        T_n = paste0(as.character(term_1), " * ", term_2)
      ) %>%
      dplyr::mutate(expr_T_n = list(parse(text = T_n))) %>%
      dplyr::ungroup()

    # Combine terms into one Taylor Series
    T_n <- call_tbl$T_n %>% paste(collapse = " + ")

    # Differentiate the last derivative we calculated
    df_np1 <- stats::D(d_list[[length(d_list)]], name = var_name)

    # Correct for interval openess if needed
    test_interval <- seq(c, d, 0.01)
    if (open_closed == "open") {
      test_interval <-
        test_interval[-which(as.character(test_interval) %in% as.character(c(c, d)))]
    }
    if (center_check) {
      test_interval <- test_interval[-which(test_interval == center)]
    }

    # Run the (n+1)th derivative's max value
    vals <-
      tibble::tibble(x = test_interval) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(f = abs(eval(df_np1, list(x = x)))) %>%
      dplyr::ungroup() %>%
      dplyr::slice_max(f) %>%
      dplyr::slice_max(x)


    # Plot the remainder function
    R_plot <-
      tibble::tibble(x = test_interval) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(R = R_nx(df_np1, x, length(d_list), a)) %>%
      dplyr::ungroup() %>%
      ggplot(aes(x, R)) + geom_line() +
      ggtitle(glue::glue("Remainder for f^{num_derivs}({var_name})"))

    # Return Taylor Series, its error, and the remainder plot
    return(list(
      T_n = T_n,
      taylor_error = (vals$f / factorial(length(d_list))) * abs(vals$x - a) ^
        (length(d_list)),
      remainder_graph = R_plot
    ))
  }

# Testing   ----

# f <- \(x) {
#   sin(x)
# }
#
# taylor <-
#   make_taylor_series(
#     f,
#     num_derivs   = 10,
#     center       = 0,
#     intr_start   = -0.5,
#     intr_end     = 0.5,
#     var_name     = "x",
#     open_closed  = "closed"
#   )
#
# taylor$taylor_error
