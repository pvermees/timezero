mdf <- function(fn){
    if (fn %in% c('He4He3_Small_Signals.mdf',
                  'He4He3_Small_Signals_Long.mdf',
                  'He4He3_Small_Signals_Very_Long.mdf')){
        out <- list(codes = c('101','102','103','104'),
                    labels = c('4He (V)','4He (cps)','3He (cps)','12C2+ (cps)'),
                    columns = c(1,2,2,2)
                    )
    } else if (fn %in% c('He4He3_Large_Signals_short.mdf',
                         'He4He3_Large_Signals.mdf',
                         'He4He3_Large_Signals_Long.mdf',
                         'He4He3_Large_Signals_Very_Long.mdf')){
        out <- list(codes = c('101','102','103'),
                    labels = c('4He (V)','3He (cps)','12C2+ (cps)'),
                    columns = c(1,2,2)
                    )
    }
    out
}
