using Distributions
using LinearAlgebra

"""
    magnitude_trigger(peak=67μV/m; dν=10.0MHz)

Create a trigger function that checks if the total magnitude of the electric field is greater than some threshold.

# Arguments
- `peak=67μV/m`: The threshold for the electric field's magnitude.
- `dν=10.0MHz`: Frequency step size.

# Returns
A function that checks if the norm of the electric field of an event is greater than the threshold.

# Description
This function generates a trigger function that evaluates whether the total magnitude of the electric field for a cosmic ray event exceeds a specified threshold. The electric field magnitude is calculated as the product of the field's polarization norm, frequency step size, and sum of electric field components.

# Usage
Used to filter events based on the strength of their electric field magnitude, allowing selection of significant cosmic ray events.
"""
function magnitude_trigger(peak=67μV / m; dν=10.0MHz)
    trig = let peak = peak, dν = dν
        event -> norm(event.pol) * dν * sum(event.Ef) >= peak
    end
    return trig
end

"""
    gaussian_trigger(peak=67μV/m, σ=45.0; dν=10.0MHz, θ0=-90.)

Creates a Gaussian beam threshold pointing down towards the nadir
    with `peak` threshold at the nadir and a width of `σ` in degrees. This uses the current estimate for a sinuous antenna array by P. Gorham

# Arguments
- `peak=67μV/m`: Peak threshold of the electric field at the nadir.
- `σ=45.0`: Width of the Gaussian beam in degrees.
- `dν=10.0MHz`: Frequency step size.
- `θ0=-90.`: Nadir angle in degrees.

# Returns
A function that checks if an event's electric field, after being scaled by a Gaussian factor based on its off-axis angle, exceeds the peak threshold.

# Description
This trigger function incorporates a Gaussian beam pattern, scaling the threshold for the electric field magnitude based on the event's off-axis angle relative to the nadir. The Gaussian factor decreases the effective threshold as the angle increases from the nadir.

# Usage
Useful in scenarios where the sensitivity of detection varies with the angle of incident cosmic rays, simulating a realistic antenna response.
"""
function gaussian_trigger(peak=67μV / m, σ=45.0; dν=10.0MHz, θ0=-90.0)
    trig = let peak = peak, σ = σ, dν = dν
        event -> (norm(event.pol) * dν * sum(event.Ef) * exp(-((θ0 - event.θ_el) * (θ0 - event.θ_el)) / (2σ * σ))) >= peak
    end
    return trig
end

"""
    rician_trigger(;snr_threshold=5.0, bw=900.0MHz, peak=calculated, σ=45.0)

Create a Rician beam trigger function based on signal-to-noise ratio.

# Parameters
- `snr_threshold=5.0`: The signal-to-noise ratio threshold for triggering.
- `bw=900.0MHz`: Bandwidth for the trigger.
- `peak=calculated`: Peak threshold calculated based on the provided parameters.
- `σ=45.0`: Width of the Gaussian beam in degrees.

# Returns
A trigger function that evaluates if the signal-to-noise ratio of an event's electric field exceeds the specified threshold.

# Description
This trigger function calculates a scaled threshold based on the event's off-axis angle and evaluates if the event's electric field, subjected to a Rician distributed noise model, surpasses this threshold. The function accounts for both primary and subsurface reflections.

# Usage
Ideal for more sophisticated detection scenarios where the signal-to-noise ratio plays a crucial role in determining event significance.
"""
function rician_trigger(; snr_threshold=5.0, bw=900.0MHz,
    peak=sqrt(85.0 / 180.0) * (1.0 / 1.9) * 67μV / m, σ=45.0)

    # calculate the equivalent threshold at this bandwidth
    # this is based on the 67μV/m at 30-300 MHz from PG
    # threshold = sqrt(bw/300.0MHz)*peak
    threshold = peak

    # the SNR at the peak that was used to calculate the 67μV/m
    peak_snr = 5.0

    # we now construct the function that implements the trigger
    trigger = let snr_threshold = snr_threshold, peak = peak, peak_snr = peak_snr, σ = σ
        function (event)

            # get the magnitude in the X-Y plane for the primary and secondary reflections
            Efmag = norm(event.Ef[1:2])
            Efmagsub = norm(event.Efsub[1:2])

            # scale the threshold up by the off-axis angle w.r.t the nadir
            scaled_threshold = threshold / exp(-(-90.0 - event.θ_el) * (-90.0 - event.θ_el) / (2σ * σ))

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



"""
    trigger_all()

Always triggers. Used for testing and debugging.
"""
function trigger_all()
    return event -> true
end

"""
    trigger_received()

Always triggers. Used for testing geometric acceptance and debugging.
"""
function trigger_received()
    return event -> !isa(event, TrialFailed)
end