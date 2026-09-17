export DLLx, DLLLx, ϵ_DLL, ϵ_DLLL, σ_DLL, umoments, ϵofx, ϵofxandt


############################################
"""
    DLLx(r; exponent_correction = 0)

Change of variables of ranges `r` to x values for a linear fit for DLL.
Optionally specify a `exponent_correction` to account for different PDF of the breakage coefficient.
"""
    function DLLx(r; exponent_correction = 0)
        return r^(2/3 + exponent_correction)  # r -> x for second-order structure function DLL
    end


"""
    DLLLx(r; exponent_correction = 0)

Change of variables of ranges `r` to x values for a linear fit for DLLL.
Optionally specify a `exponent_correction` to account for different PDF of the breakage coefficient.
"""
    function DLLLx(r; exponent_correction = 0)
        return r^(1 + exponent_correction)  # r -> x for third-order structure function DLLL
    end


"""
    ϵ_DLL(A; C₂ = 2.0)

Calculate dissipation rate ϵ from the slope `A` of a linear fit to DLL.
"""
    function ϵ_DLL(A; C₂ = 2.0)
        return sign(A) * abs(A / C₂)^(3/2) # Estimate from DLL
    end


"""
    ϵ_DLLL(A; C₂ = 2.0)

Calculate dissipation rate ϵ from the slope `A` of a linear fit to DLLL.
"""
    function ϵ_DLLL(A)
        return -5/4 * A # Estimate from DLLL
    end


"""
    σ_DLL(B)

Calculate the noise variance from the intercept of the fit to DLL.
"""
    function σ_DLL(B)
        return B / 2
    end


"""
    TimeDuration(t)

Calculate the duration of the time series `t`.
"""
    function TimeDuration(t::AbstractVector{<:Real})
        return t[end] - t[1]
    end

"""
    TimeDuration(t)

Calculate the duration of the time series `t` in seconds.
"""
    function TimeDuration(t::AbstractVector{<:TimeType})
        return Dates.value(t[end] - t[1])/1000
    end


############################################
### Calculating the Velocity Moments ###

"""
    umoments(ur, Δr, x, rmax; min_points = 4, calculate_sf3 = true)

Calculate the second moments of longitudinal velocity `ur` with separation `Δr` and coordinate
`x` using all possible separations in the region `rmax` around `x`.

`min_points` is the minimum number of values of `ur` to to use to calculate the velocity moments (default is 4, resulting in 3 separations).
Optionally calculate the third moments of `ur` by setting `calculate_sf3 = true` (default).
"""
function umoments(ur, Δr, x, rmax; min_points = 4, calculate_sf3 = true)
    
        ## Determine the separations and the maximum size of the window ##
        Δx = mean(diff(x))       # x separation
        ΔN = round(Int, rmax/Δr, RoundNearestTiesUp) # The number of points in a window
        #####

        ### Pad the arrays to be able to get the edges ###
        ####################################
        ## The arrays are padded such that the estimates near the edges have `min_points` of velocity measurements ##
        # the third term is to ensure the estimates are within the domain by removing padding that would allow an outside estimate 
        pad_size = ΔN - (min_points - 1) - ceil.(Int, ΔN./2 .- 2) # the number of points to pad with

        padding = 1:pad_size      # Initialize the padding vector as increasing to be able to extend the coordiantes

        ## a continuation of the coordinate z on either side of the vector ##
        pad_x_start = x[1] .- Δx*reverse(padding)
        pad_x_end = x[end] .+ Δx*padding

        ## Need to pad `ur` in 2 dimensions
        Nt = size(ur, 2) # the number of time steps
        ur_padding = fill!(Array{eltype(ur)}(undef, pad_size, Nt), NaN)

        ## Pad the arrays ##
        x_pad = cat(pad_x_start, x, pad_x_end, dims=1)
        ur_pad = cat(ur_padding, ur, ur_padding, dims=1)

        ## Make a separation vector for algebra. The actual values don't matter, only their difference ##
        r = collect(Δr.*(1:size(x, 1)))
        r_pad = cat(NaN .* padding, r, NaN .* padding, dims=1)

        #####################################

        N = size(x_pad, 1) - ΔN # The number of possible estimates

        ## The coordinate vector will be the same size as the input
        x_assigned = fill!(similar(x), NaN)

        ML2 = Array{Array{eltype(ur)}}(undef, N)
        if calculate_sf3
            ML3 = Array{Array{eltype(ur)}}(undef, N)
        end

        ML2_nomean = Array{Array{eltype(ur)}}(undef, N)
        if calculate_sf3
            ML3_nomean = Array{Array{eltype(ur)}}(undef, N)
        end


        @threads for r_idx in 1:N

            ## The indices of the box around the position `r_idx` ##
            box_idx = r_idx:(r_idx + ΔN)

            ## Give the point a coordinate position ##
            x_assigned[r_idx] = mean(x_pad[box_idx])

            ## Find all the possible separations ##
            Δr_mat = r_pad[box_idx] .- r_pad[box_idx]'

            pos_idx = findall(x-> x>0, Δr_mat) # take only the positive separations

            ## Create a velocity matrix with an extra empty dimension ##
            w_zt = reshape(ur_pad[box_idx,:], ΔN+1, 1, Nt)

            ## Find all possible velocity differences ##
            Δw_zt = w_zt .- permutedims(w_zt, [2,1,3]) # The extra empty dimension allows for a higher dimensional transpose to make this line efficient


            ## Square the velocity differences, then take their mean and add them to the second column of the array ##
            pos_r = Δr_mat[pos_idx]

            ML2_n = Δw_zt.^2
            ML2_nomean[r_idx] = vec(ML2_n[pos_idx,:]) # Add the non-meaned data (for QA) to a different data structure
            ML2[r_idx] = hcat(pos_r, nanmean(ML2_n; dims=3)[pos_idx])

            if calculate_sf3
            ## And cube them ##
                ML3_n = Δw_zt.^3
                ML3_nomean[r_idx] = vec(ML3_n[pos_idx,:]) # Add the non-meaned data (for QA) to a different data structure
                ML3[r_idx] = hcat(pos_r, nanmean(ML3_n, dims=3)[pos_idx])
            end

        end

        # Put everything into one named tuple
        if calculate_sf3
            return (; x = x_assigned, ML2, ML2_nomean, ML3, ML3_nomean, true_rmax = Δr * ΔN)
        else
            return (; x = x_assigned, ML2, ML2_nomean, true_rmax = Δr * ΔN)
        end
    end


############################################

### Calculating TKE Dissipation Rate ###

"""
    ϵofx(umom; fittype=:OLS, γ = 0.32, minp2fit = 3, exponent_correction = 0, calculate_sf3 = true)

Calculate TKE dissipation rate ϵ from the velocity moments returned by `umoments`, and the doppler noise variance.

Calculate from both the second and third order moments.

...
# Arguments
- `fittype`: The method used for fitting the line (:OLS or :GLS).
- `γ`: The confidence level (i.e. γ-th percentile).
- `minp2fit`: The minimum number of points to fit to.
- `exponent_correction`: optional correction factor for the exponent.
- `calculate_sf3`: Whether to calculate ϵ from the third-order moments (default is true).
...
"""
    function ϵofx(umom; fittype=:OLS, γ = 0.32, minp2fit = 3, exponent_correction = 0, calculate_sf3 = true)

        N = size(umom.ML2, 1) # Number of points to estimate ϵ at

        N_coord = size(umom.x, 1) # Size of the input/output array

        ## Set the least squares fit parameters ##
        fitfunc = getfield(FourierFuncs, fittype) # get the least squares fit function (OLS or GLS)
        fmodel(x, p) = p[1] * x .+ p[2]
        fmodel_power(x, p) = p[1] * x .^ p[2]

        ## Initialize the arrays to fill ##
        data_type = eltype(umom.ML2[1]) #Float64 # get the data type from ML2
        data_fill = NaN

        sf2_ϵ = fill!(Array{data_type}(undef, N_coord), data_fill)
        sf2_ϵ_uncertainty = fill!(Array{data_type}(undef, N_coord), data_fill)
         sf2n = fill!(Array{data_type}(undef, N_coord), data_fill) # For the best fit power QA

         vari = fill!(Array{data_type}(undef, N_coord), data_fill)
        vari_uncertainty = fill!(Array{data_type}(undef, N_coord), data_fill)

        if calculate_sf3
            sf3_ϵ = fill!(Array{data_type}(undef, N_coord), data_fill)
            sf3_ϵ_uncertainty = fill!(Array{data_type}(undef, N_coord), data_fill)
             sf3n = fill!(Array{data_type}(undef, N_coord), data_fill) # For the best fit power QA
        end


        ## Loop through each x coordinate in `umom`
        @threads for x_idx in 1:N

            ### Extract the variables ###
            ## Second-order ##
            ML22fit = umom.ML2[x_idx][:,2]          # The second moments
            sf2x2fit = umom.ML2[x_idx][:,1]         # The separation axis to fit

            ## Third-order ##
            if calculate_sf3
                ML32fit = umom.ML3[x_idx][:,2]          # The Third moments
                sf3x2fit = umom.ML3[x_idx][:,1]         # The separation axis to fit
            end


            ### Find the non-NaN values to fit to ###
            # The NaNs will be consistent between ML2 and ML3
            p2fit = .!isnan.(ML22fit)            # Find the non-NaN values to fit to

            if sum(p2fit) >= minp2fit # make sure there is a minimum number of point used to fit

                ### Least squares fit the velocity moments as a function of separation ###
                LSfit_DLL = fitfunc(DLLx.(sf2x2fit[p2fit]; exponent_correction), ML22fit[p2fit]; fitmodel = fmodel) # fit


                ## Calculate the confidence intervals ##
                # SF2 #
                sf2_Δ = linear_uncertainty(DLLx.(sf2x2fit[p2fit]; exponent_correction), ML22fit[p2fit], fmodel(DLLx.(sf2x2fit[p2fit]; exponent_correction), LSfit_DLL.param), γ)


                ## Calculate ϵ from the slopes of the fits and the confidence intervals ##
                sf2_ϵ_measurement = ϵ_DLL(LSfit_DLL.param[1] ± sf2_Δ.beta_e)
                    sf2_ϵ[x_idx] = Measurements.value(sf2_ϵ_measurement)
                    sf2_ϵ_uncertainty[x_idx] = Measurements.uncertainty(sf2_ϵ_measurement)

                vari_measurement = σ_DLL(LSfit_DLL.param[2] ± sf2_Δ.alpha_e) # The doppler noise from the intercept
                    vari[x_idx] = Measurements.value(vari_measurement)
                    vari_uncertainty[x_idx] = Measurements.uncertainty(vari_measurement)


                ### Calculate the best fit power for QA ###
                DLL_power = fitfunc(sf2x2fit[p2fit], ML22fit[p2fit] .- LSfit_DLL.param[2]; fitmodel = fmodel_power)
                sf2n[x_idx] = DLL_power.param[2]


                if calculate_sf3 # optionally calculate ϵ from the third-order moments
                    LSfit_DLLL = fitfunc(DLLLx.(sf3x2fit[p2fit]), ML32fit[p2fit]; fitmodel = fmodel) # fit
                    sf3_Δ = linear_uncertainty(DLLLx.(sf3x2fit[p2fit]), ML32fit[p2fit], fmodel(DLLLx.(sf3x2fit[p2fit]), LSfit_DLLL.param), γ)
                    sf3_ϵ_measuement = ϵ_DLLL(LSfit_DLLL.param[1] ± sf3_Δ.beta_e)
                        sf3_ϵ[x_idx] = Measurements.value(sf3_ϵ_measuement)
                        sf3_ϵ_uncertainty[x_idx] = Measurements.uncertainty(sf3_ϵ_measuement)
                    DLLL_power = fitfunc(sf3x2fit[p2fit], ML32fit[p2fit] .- LSfit_DLLL.param[2]; fitmodel = fmodel_power)
                    sf3n[x_idx] = DLLL_power.param[2]
                end
            end
        end

        if calculate_sf3
            return (; sf2_ϵ, sf2_ϵ_uncertainty, sf3_ϵ, sf3_ϵ_uncertainty, vari, vari_uncertainty, sf2n, sf3n)
        else
            return (; sf2_ϵ, sf2_ϵ_uncertainty, vari, vari_uncertainty, sf2n)
        end
    end



"""
    ϵofxandt(ur, x, t, Δr, fs, Δt, rmax; xmin = minimum(x), xmax = maximum(x), min_points = 4, fittype = :OLS, γ = 0.32, minp2fit = 3, calculate_sf3 = true)

Calculate a TKE dissipation rate time series using velocities `ur` with separation `Δr`, spatial coordinate `x`, times `t`, and downsample time resolution `Δt`.

...
# Arguments
- `ur`: The matrix of velocities.
- `x`: The spatial coordinates of `ur`.
- `t`: The times of `ur`.
- `Δr`: The separation distance between velocities.
- `fs`: the sample frequency of `t`.
- `Δt`: the time interval to average velocity moments over. This becomes the resulting time resolution.
- `rmax`: the maximum range separation between velocities to consider.
- `xmin`: The minimum value in `x` to include in the calculations. Can either be a single value or a vector of `size(t)`.
- `xmax`: The maximum value in `x` to include in the calculations. Can either be a single value or a vector of `size(t)`.
- `min_points`: The minimum number of values of `ur` to use to calculate ϵ (default is 4, resulting in 3 separations).
- `fittype`: The method used for fitting the line (:OLS or :GLS).
- `γ`: The confidence level (i.e. γ-th percentile).
- `minp2fit`: The minimum number of points to fit to.
- `exponent_correction`: optional correction factor for the exponent.
- `calculate_sf3`: Whether to calculate ϵ from the third-order moments (default is true).
...
"""
    function ϵofxandt(ur, x, t, Δr, fs, Δt, rmax; xmin = minimum(x), xmax = maximum(x), min_points = 4, fittype = :OLS, γ = 0.32, minp2fit = 3, exponent_correction = 0, calculate_sf3 = true)

        ## Make a vector of indices to have the correct number of segments in the time series ##
        time_idx = 1:floor(Int, TimeDuration(t)/Δt)

        ## Make the spatial coordinate vector ##
        ΔN = round(Int, rmax/Δr, RoundNearestTiesUp) # The number of points in a window
        Δx = mean(diff(x))

        if iseven(ΔN)
            x_coords = x
        else
            x_coords = x .+ Δx/2
        end


        ## Find which indices of the data to use depending on the depth ##
        x_mat = repeat(x, 1, size(t, 1))
        x_idx_t = xlims2idx(x_mat, xmin, xmax)


        ## Initialize some arrays to fill ##
               Nt = size(time_idx, 1) # Number of times
               Nz = size(x, 1)        # Number of positions in the `x` coordinate
        data_type = eltype(ur) # get the data type from ur
        data_fill = NaN

              times = Array{eltype(t)}(undef, Nt)
        umomofxandt = Array{NamedTuple}(undef, Nt)
                sf2 = fill!(Array{data_type}(undef, Nz, Nt), data_fill)
               sf2σ = fill!(Array{data_type}(undef, Nz, Nt), data_fill)
               vari = fill!(Array{data_type}(undef, Nz, Nt), data_fill)
              variσ = fill!(Array{data_type}(undef, Nz, Nt), data_fill)
               sf2n = fill!(Array{data_type}(undef, Nz, Nt), data_fill)

        if calculate_sf3
             sf3 = fill!(Array{data_type}(undef, Nz, Nt), data_fill)
            sf3σ = fill!(Array{data_type}(undef, Nz, Nt), data_fill)
            sf3n = fill!(Array{data_type}(undef, Nz, Nt), data_fill)
        end

        ## Now loop through each time ##
        @threads for tt in time_idx

            # Make the vector of indices for the time interval
            tt_idx = Int.(fs*Δt*(tt-1)+1:fs*Δt*tt)

            # Find the center of the interval for assigning a time
            time_cen = round(Int, median(tt_idx), RoundNearestTiesUp) # round so that it is the closest time, need to specify the rounding mode for it to work

            ## Calculate moments ##
            umom = umoments(ur[x_idx_t[time_cen], tt_idx], Δr, x_mat[x_idx_t[time_cen], time_cen], rmax; min_points)

            ## Calculate ϵ from moments ##
            ϵ = ϵofx(umom; fittype, γ, minp2fit, exponent_correction, calculate_sf3)

            ## Store the times ##
            times[tt] = t[time_cen]

            ## Keep the moments for reference ##
            umomofxandt[tt] = umom

            ## Store the ϵ values ##
            sf2[x_idx_t[time_cen], tt] = ϵ.sf2_ϵ # value
            sf2σ[x_idx_t[time_cen], tt] = ϵ.sf2_ϵ_uncertainty # uncertainty

            ## Store the doppler noise variance from SF2 ##
            vari[x_idx_t[time_cen], tt] = ϵ.vari # value
            variσ[x_idx_t[time_cen], tt] = ϵ.vari_uncertainty # uncertainty

            ## Store the power metrics ##
                sf2n[x_idx_t[time_cen], tt] = ϵ.sf2n
    
            if calculate_sf3
                sf3[x_idx_t[time_cen], tt] = ϵ.sf3_ϵ
                sf3σ[x_idx_t[time_cen], tt] = ϵ.sf3_ϵ_uncertainty # uncertainty
                sf3n[x_idx_t[time_cen], tt] = ϵ.sf3n
            end
        end

        ## display the actual maximum separation ##
        true_rmax = umomofxandt[1].true_rmax
        show("Actual maximum separation = $(true_rmax)")
        ####

        if calculate_sf3
            return (; x=x_coords, t=times, sf2, sf2σ, vari, variσ, sf2n, sf3, sf3σ, sf3n, umom = umomofxandt, true_rmax)
        else
            return (; x=x_coords, t=times, sf2, sf2σ, vari, variσ, sf2n, umom = umomofxandt, true_rmax)
        end
    end