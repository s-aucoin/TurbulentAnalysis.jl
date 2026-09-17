export FailsQC, QCϵts


############################################
"""
    FailsQC(q, n)
Determine if a point with QC flag `q` fails quality control condition number `n`.
"""
    function FailsQC(q, n)
        return ((q >> n) & 1) == 1
    end


"""
    QCϵts(ϵ_ts, max_relative_uncertainty, max_relative_power_diff; exponent_correction = 0)

Replace values of ϵ in `ϵ_ts` that exceed the `max_relative_uncertainty` or `max_relative_power_diff` with `fillvalue`.
Optionally specify a `exponent_correction` to account for different PDF of the breakage coefficient.
"""
    function QCϵts(ϵ_ts, max_relative_uncertainty, max_relative_power_diff; exponent_correction = 0)

        ## Check if there are any NaNs (from the position being out of bounds when calculated) ##
        nan_sf2 = isnan.(ϵ_ts.sf2)
        nan_vari = isnan.(ϵ_ts.vari)
        if haskey(ϵ_ts, :sf3)
            nan_sf3 = isnan.(ϵ_ts.sf3)
        end

        ## find the points with unphysical estimates (ϵ < 0) ##
        negative_sf2 = ϵ_ts.sf2 .< 0
        negative_vari = ϵ_ts.vari .< 0
        if haskey(ϵ_ts, :sf3)
            negative_sf3 = ϵ_ts.sf3 .< 0
        end

        ## find the points with unacceptable uncertainties ##
        bad_data_sf2 = (ϵ_ts.sf2σ ./ ϵ_ts.sf2) .> max_relative_uncertainty
        bad_data_vari = (ϵ_ts.variσ ./ ϵ_ts.vari) .> max_relative_uncertainty
        if haskey(ϵ_ts, :sf3)
            bad_data_sf3 = (ϵ_ts.sf3σ ./ ϵ_ts.sf3) .> max_relative_uncertainty
        end

        ## Find the points with unacceptable best fit powers##
        sf2_expected_power = 2/3
        sf2_max_power_diff = (sf2_expected_power + exponent_correction) * max_relative_power_diff

        bad_power_sf2 = (ϵ_ts.sf2n) .< (sf2_expected_power - sf2_max_power_diff) .|| 
                        (ϵ_ts.sf2n) .> (sf2_expected_power + sf2_max_power_diff)

        if haskey(ϵ_ts, :sf3)
            sf3_expected_power = 1
            sf3_max_power_diff = (sf3_expected_power + exponent_correction) * max_relative_power_diff

            bad_power_sf3 = (ϵ_ts.sf3n) .< (sf3_expected_power - sf3_max_power_diff) .|| 
                            (ϵ_ts.sf3n) .> (sf3_expected_power + sf3_max_power_diff)
        end

        ## assign each quality control condition a bit ##
        nan_flag = Int8(2^0)
        negative_flag = Int8(2^1)
        uncertainty_flag = Int8(2^2)
        power_flag = Int8(2^3)

        qsf2 = (nan_flag*nan_sf2) + (negative_flag*negative_sf2) + (uncertainty_flag*bad_data_sf2) + (power_flag*bad_power_sf2)
        qvari = (nan_flag*nan_vari) + (negative_flag*negative_vari) + (uncertainty_flag*bad_data_vari)
        if haskey(ϵ_ts, :sf3)
            qsf3 = (nan_flag*nan_sf3) + (negative_flag*negative_sf3) + (uncertainty_flag*bad_data_sf3) + (power_flag*bad_power_sf3)
        end

        if haskey(ϵ_ts, :sf3)
            return (; qsf2, qvari, qsf3)
        else
            return (; qsf2, qvari)
        end
    end