using DelimitedFiles
using Unitful: Ω

"""
Construct an antenna-factor based antenna simulation.
"""
function create_antenna(fname, skyfrac;
                        altitude=20.0km,
                        ν_min=150MHz, ν_max=800MHz,
                        σθ=39.1, σϕ=28.4, θ0=-40.0,
                        Nant=8, Ntrig=3, Tlna=60K, Tmoon=85K, SNR=4.0, )

    # current LPDA design is 150-600 MHz
    # ANITA horns are 100-800 MHz
    ν_ant = ν_min:10MHz:ν_max

    # the antenna-receiver matching
    Γ = 0.9

    # the fraction of sky and moon that the LPDA's see
    moonfrac = 1.0 - skyfrac

    # get the noise temperature of the sky
    Tsky = sky_temperature(ν_ant)

    # construct the array of angles to calculate the antenna temperature
    θb = 0.0:0.5:180.0

    # construct the beam weightings at each angle
    # we analytically do the ϕ integral so it's (2π)*(1.0 / 4π) = 0.5
    Rb = @. 0.5 * exp(-((90.0 - θ0) - θb)*((90.0 - θ0) - θb) / (2σθ*σθ))

    # and now calculate the integral over θ and frequency
    # and then sum over angles
    # this is now the antenna temperature, θb[2] - θb[1] is dθ
    Tbeam = @.( Rb *sind(θb) * deg2rad(θb[2] - θb[1])) * Tsky'

    # calculate the horizon angle at this altitude
    θh = 90.0 + rad2deg(-horizon_angle(altitude))

    # and replace the points below the horizon with the moon temperature
    Tbeam[θb .> θh, :] .= (Rb .* Tmoon .* sind.(θb) * deg2rad(θb[2] - θb[1]))[θb .> θh]

    # and sum it over the angles
    Tant = sum(Tbeam, dims=1)

    # construct the system temperature
    Tsys = Tlna .+ Γ*(Tsky.*skyfrac .+ Tmoon*moonfrac)
    #Tsys = Tlna .+ Γ*Tant

    # the current LPDA design is nominally a 120 ohm design
    # but the antenna factors are referenced to 50 ohm so
    # calculate the noise at 50 ohm directly
    Z = 50Ω

    # load the antenna factor from files
    AFs = readdlm("$(@__DIR__)/../data/$(fname).txt", skipstart=1)

    # create an interpolator and interpolate it onto our frequencies
    AF = LinearInterpolation(AFs[:, 1] ./ 1e6, AFs[:, 2])(ν_ant ./ MHz) ./ m

    # construct the noise voltage spectrum - this is units of Volts
    Vn = sqrt.(2.0 * Ntrig * k_b .* Tsys .* Z * (ν_ant[2] - ν_ant[1])) .|> μV

    # we also need the total summed noise voltage
    Vntot::typeof(1.0μV) = sqrt(sum(Vn.^2))

    # construct the azimuthal angles given the number of antennas
    ϕ = range(0, 2π, length=Nant+1)[1:Nant]

    # construct the boresight vectors for each antenna
    boresight = spherical_to_cartesian.(π/2.0 .+ deg2rad(-θ0)*ones(Nant), ϕ, 1)

    # construct the perpendicular (i.e. H-pol) antenna axis
    Haxis = [boresight[i] × SA[0, 0, 1] for i=1:Nant]
    Haxis ./= norm.(Haxis) # make sure it is normalized

    # and now cross the Haxis with the boresight to get the Vaxis
    Vaxis = Haxis .× boresight
    Vaxis ./= norm.(Vaxis) # make sure it is normalized

    # we now construct the function that evaluates the trigger
    # use a let-block to improve performance on earlier Julia versions
    trigger = let ν_ant=ν_ant, Vntot=Vntot, AF=AF, Nant=Nant,
        Ntrig=Ntrig, Haxis=Haxis, Vaxis=Vaxis, boresight=boresight
        function (event::AbstractSignal)

            # make sure our simulated LPDA frequencies are
            # the same as the electric field frequencies
            @assert (ν_ant[2] - ν_ant[1]) == event.dν

            # interpolate the electric field onto the antennas frequencies
            # VmMHz = (1V/m/MHz) .* LinearInterpolation((event.ν_min:event.dν:event.ν_max) ./ MHz,
            #                                           event.Ef ./ (V/m/MHz))(ν_ant ./ MHz)
            VmMHz = (1μV/m/MHz) .* extrapolate(interpolate(((event.ν_min:event.dν:event.ν_max) ./ MHz,),
                                              event.Ef ./ (μV/m/MHz),
                                                           Gridded(Linear())), 0.0)(ν_ant ./ MHz)

            # and calculate the electric field spectrum in V/m
            Ef = norm(event.pol) * VmMHz .* (ν_ant[2] - ν_ant[1]) # now in Volts/m

            # construct the vector to the event from SC
            view = spherical_to_cartesian.(π/2.0 + deg2rad(-event.θ_el),
                                           deg2rad(event.ϕ), 1)

            # this array stores the relative antenna response factors due
            # to the off-axis beam response
            scale = ones(Nant)

            # loop over every antenna to calculate θ, ϕ and then the off-axis weighting
            for i = 1:Nant

                # project the view vector onto the plane defined by Haxis
                # i.e. the vertical plane centered at the boresight
                Vproj = view - (view ⋅ Haxis[i])*Haxis[i]

                # the angle between `Vproj` and boresight is the elevation angle
                θev = acos( (Vproj / norm(Vproj)) ⋅ boresight[i] ) |> rad2deg

                # project the view vector onto the plane defined by Vaxis
                # i.e. the horizontal plane centered at the boresight
                Hproj = view - (view ⋅ Vaxis[i])*Vaxis[i]

                # the angle between `Hproj` and boresight is the azimuth angle
                ϕev = acos( (Hproj / norm(Hproj)) ⋅ boresight[i] ) |> rad2deg

                # apply the off-axis antenna response at these angles
                scale[i] *= 0.074 + (1. - 0.074)exp(-θev*θev / (2σθ*σθ))
                scale[i] *= 0.062 + (1. - 0.062)exp(-ϕev*ϕev / (2σϕ*σϕ))

            end

            # we now pick the best Ntrig antennas to calculate the trigger with
            Nsum = sum(sort!(scale, rev=true)[1:Ntrig])

            # apply the antenna factors to the electric field
            # to convert it to Volts, and scale it by the number of antennas
            # double the field because we have two polarizations per antenna
            Vs = (2.0 * Nsum) * (Ef ./ AF)

            # and lastly sum up the electric field
            Vf = sum(Vs)

            # do the same process above for the reflection from the bottom of the layer
            if typeof(event) == Reflected

                # and calculate the electric field spectrum for the bottom reflection in V/m
                Efsub = norm(event.polsub) * VmMHz .* (ν_ant[2] - ν_ant[1]) # now in Volts/m

                # and calculate the corresponding voltage at the receiver
                Vssub = (2. * Nsum) * (Efsub ./ AF)

                # and sum it up for the trigger
                Vfsub = sum(Vssub)
            else
                Vfsub = 0.0μV
            end

            # realize a random signal+noise spectrum
            # this is in volts already, not V/MHz
            # throw for two different events - surface and subsurface
            # Vf = rician(Vs ./ μV, Vn ./ μV) .* μV
            # Vfsub = rician(Vssub ./ μV, Vn ./ μV) .* μV

            # throw for a Rician with signal and noise
            # println(sum(Vs))
            # println(Vntot)
            # Vf = rician( sum(Vs) / μV,  Vntot / μV)μV
            # println(Vf)
            # Vfsub = 0.0μV

            # println(sum(Vs) / Vntot)

            # println(sum(Vs))
            # println(sum(Vssub))
            # println(Vntot)

            # check that the total integrated voltage is greater
            # than the SNR times the summed noise voltage
            # we only need to check the greater of the top
            # and bottom surface reflections
            # return max(sqrt(sum(Vf.^2)), sqrt(sum(Vfsub.^2))) > (SNR*Vntot)
            return max(Vf, Vfsub) > (SNR*Vntot)
            # return sum(Vs) > (SNR * Vntot)
            #return max(sum(Vs), sum(Vssub)) > (SNR*Vntot)
            # return max(Vf, Vfsub) > (SNR*Vntot)

        end # end trigger function
    end # end let-block

    return trigger
end # end create_antenna

# create a model for LPDA's and for ANITA horns
LPDA(;kwargs...) = create_antenna("dualLPDA_antennaFacsV3", 0.03; kwargs...)
ANITA(;kwargs...) = create_antenna("ANITA_antennaFacs", 0.03; ν_min = 150MHz, ν_max=800MHz, kwargs...)
