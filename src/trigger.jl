using Distributions
using LinearAlgebra

"""
    magnitude_trigger(threshold)

Implement a simple trigger that checks the total magnitude
of the electric field is greater than some threshold.

This returns a *function* - that function checks if the norm of
the electric field of that event is greater than some pre-defined
threshold.
"""
function magnitude_trigger(peak=67μV/m; dν=10.0MHz)
    trig = let peak=peak, dν = dν
        event -> norm(event.pol)*dν*sum(event.Ef) >= peak
    end
    return trig
end

"""
    gaussian_beam(peak, σ)

Implement a Gaussian beam threshold pointing down towards the nadir
with `peak` threshold at the nadir and a width of `σ` in degrees.

This uses the current estimate for a sinuous antenna array by P. Gorham
"""
function gaussian_trigger(peak=67μV/m, σ = 45.0; dν=10.0MHz, θ0 = -90.)
    trig = let peak=peak, σ=σ, dν = dν
        event -> (norm(event.pol)*dν*sum(event.Ef)*exp( -((θ0 - event.θ_el)*(θ0 - event.θ_el))  / (2σ*σ))) >= peak
    end
    return trig
end

"""
    rician_beam(threshold)

Implement a Gaussian beam threshold pointing down towards the nadir.

This uses the current estimate for a sinuous antenna array by P. Gorham
"""
function rician_trigger(;snr_threshold=5.0, bw=900.0MHz,
                        peak=sqrt(85.0 / 180.0)*(1.0/1.9)*67μV/m, σ = 45.0)

    # calculate the equivalent threshold at this bandwidth
    # this is based on the 67μV/m at 30-300 MHz from PG
    # threshold = sqrt(bw/300.0MHz)*peak
    threshold = peak

    # the SNR at the peak that was used to calculate the 67μV/m
    peak_snr = 5.0

    # we now construct the function that implements the trigger
    trigger = let snr_threshold=snr_threshold, peak=peak, peak_snr=peak_snr, σ=σ
        function (event)

            # get the magnitude in the X-Y plane for the primary and secondary reflections
            Efmag = norm(event.Ef[1:2])
            Efmagsub = norm(event.Efsub[1:2])

            # scale the threshold up by the off-axis angle w.r.t the nadir
            scaled_threshold = threshold / exp( -(-90.0 - event.θ_el)*(-90.0 - event.θ_el) / (2σ*σ) )

            # throw for a Rician SNR threshold - centered at SNR of 5.0
            # since this is the SNR PG used to calculate the 67μV/m
            # we do this for the primary event and the subsurface event
            snr = rician(peak_snr * Efmag / scaled_threshold, 1)
            snrsub = rician(peak_snr * Efmagsub / scaled_threshold, 1)

            # if this is a direct event, ignore the calculation of the subsurface trig
            if event.type == DirectEvent
                return snr > snr_threshold
            else
                # allow for either trigger
                return (snr > snr_threshold) || (snrsub > snr_threshold)
            end

            # we can never get here
        end
    end

    return trigger
end
