export largest_eddies


############################################
"""
    largest_eddies(u, x, t, Δr, δt, fs; xmin = minimum(x), xmax = maximum(x), x0_c = nothing, SF = nothing, logmean=false)

Estimate the properties of the largest eddies of `u`.

...
# Arguments
- `u`: longitudinal velocities along the line defined by `x`. `t` is the second dimension.
- `x`: The spatial coordinates of `u`.
- `t`: The times of `u`.
- `Δr`: The separation distance between velocities.
- `δt`: the time interval to average over.
- `fs`: the sample frequency of `t` (units of 1/δt).
- `xmin`: The minimum value in `x` to include in the calculations. Can either be a single value or a vector of `size(t)`.
- `xmax`: The maximum value in `x` to include in the calculations. Can either be a single value or a vector of `size(t)`.
- `x0_c`: optional reference position to calculate covariance from.
- `SF`: optional structure function ϵ object to average to compare.
- `logmean`: Whether to take the logarithmic mean of the `SF` ϵ estimates.
...
"""
    function largest_eddies(u, x, t, Δr, δt, fs; xmin = minimum(x), xmax = maximum(x), x0_c = nothing, SF = nothing, logmean=false)

        δt_idx = Int(fs*δt / 2)  # convert the size of the window to index

        δx = xlims2idx(repeat(x, 1, size(u)[2]), xmin, xmax)  # indices of the depth interval to use

        t_idx = 1 + δt_idx:size(u)[2] - δt_idx                # The indices of the times to calculate at


        times = t[t_idx]

        # Initialize some arrays #
            ell = Array{Float64}(undef, length(t_idx))
          u_rms = Array{Float64}(undef, length(t_idx))
        ϵ_rough = Array{Float64}(undef, length(t_idx))
          ϵ_sf2 = Array{Float64}(undef, length(t_idx))
          ϵ_sf3 = Array{Float64}(undef, length(t_idx))

        # Loop through each time #
        @threads for (iter, t_cen) in collect(enumerate(t_idx))

            t_sub = t_cen - δt_idx:t_cen + δt_idx # the indices of the time window

            ## Remove the mean profile ##
            u_sub = u[:,t_sub] .- mean(u[:,t_sub], dims=2)

        # If a single reference x0 is chosen #
        ########################################
            if x0_c != nothing
                # Calcuate the covariance of each depth with the reference depths
                covmat = Array{Float64}(undef, length(δx[iter]), length(t_sub))
                @threads for x in δx[iter]                             # loop through each depth
                    covmat[x,:] = u_sub[x,:] .* u_sub[x0_c,:] # multiply element-wise to get u(x,t)u(x0,t)
                end
                covvec = vec(mean(covmat, dims=2))            # average in time
        ########################################

        # If not, average over each possible x0 #
            else
                # Calcuate the covariance of each depth with the reference depths
                covmat = Array{Float64}(undef, length(δx[iter]), length(t_sub), length(δx[iter]))
                @threads for (x, x0) in collect(Iterators.product(δx[iter], δx[iter])) # loop through each depth
                    covmat[x,:,x0] = u_sub[x,:] .* u_sub[x0,:]                # multiply element-wise to get u(x,t)u(x0,t)
                end
                covvec = vec(mean(covmat, dims=(2,3)))                        # average in time and reference depth to get a <u(x)u(x0)>
            end

        ########################################

            # Calculate the rms velocity of the time window
            u_rms_n = sqrt(mean(u_sub[δx[iter],:].^2))
            u_rms[iter] = u_rms_n


            # Convert the covariance to correlation and integrate to get the lengthscale
            ell_n = int_def_trap(covvec, Δr) / (u_rms_n.^2)
            ell[iter] = ell_n

            # estimate the TKE dissipation rate
            ϵ_rough[iter] = u_rms_n^3/ell_n

            # Optionally find the equivalent from the structure function #
            if SF != nothing
                # Find the equivalent time interval #
                SFtimes = SF.t
                 (val, stidx) = findmin(abs.((SFtimes .- SFtimes[1]) .- t_sub[1]/fs))
                (val, endidx) = findmin(abs.((SFtimes .- SFtimes[1]) .- t_sub[end]/fs))

                if logmean
                    ϵ_sf2[iter] = 10 .^ nanmean(log10.(SF.sf2[:,stidx:endidx]))
                    ϵ_sf3[iter] = 10 .^ nanmean(log10.(SF.sf3[:,stidx:endidx]))
                else
                    ϵ_sf2[iter] = nanmean(SF.sf2[:,stidx:endidx])
                    ϵ_sf3[iter] = nanmean(SF.sf3[:,stidx:endidx])
                end
            end

        end

        # remove the negative estimates #
        ϵ_rough[ϵ_rough .< 0] .= NaN

        return (; times, ell, u_rms, ϵ_rough, ϵ_sf2, ϵ_sf3)

    end