using LinearAlgebra
using Statistics
using JLD2
using MATLAB
using ProgressMeter
@load joinpath(dirname(@__FILE__),"orbit_data_zac.jld2") τ_hist B_hist_b J

function run_single_sim(τ_hist,B_hist_b,J,W_telescope)
    """Run a single simulation with specified telescope sensing error"""
    dt = 1.0 #seconds
    t = Array(1:600)

    J_boom = Diagonal(0.25.*[.1^2; .2^2; .3^2])
    rad2arcsec = 3600*180/pi

    Nm2μNm = 1e6
    τ = Nm2μNm.*hcat(τ_hist...)

    #Double integrator dynamics
    #Units are arcsec, seconds, torque in μNm
    #State is [θ; θ̇]
    A = [zeros(3,3) I; zeros(3,6)]
    B = Array([zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J)]) #control input jacobian
    C = Array(Diagonal(ones(6)))
    G = [zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J)] #disturbance input jacobian

    #Convert to discrete time
    H = exp(dt.*[A B G; zeros(6,12)])
    A = H[1:6,1:6]
    B = H[1:6,7:9]
    G = H[1:6,10:12];

    #Kalman Filter Design
    V = [0.00001*I zeros(6,3); zeros(3,6) 0.0001*I] # TODO: tune this


    # W_telescope = 0.001^2 #arcsec^2 1-sigma at 1 Hz
    #W_gyro = (0.06*60)^2 #(arcsec/sec)^2 1-sigma at 1 Hz for Epson MEMS IMU (ARW = 0.06 deg/sqrt(hr))
    W_gyro = (0.0035*60)^2 #(arcsec/sec)^2 1-sigma at 1 Hz for Honeywell GG1320 laser gyro (ARW = 0.0035 deg/sqrt(hr))
    W = Array(Diagonal([(W_telescope/3.0)*ones(3); (W_gyro/3.0)*ones(3)])) #

    #Double integrator dynamics + bias torque
    #Units are arcsec, seconds, torque in μNm
    #State is [θ; θ̇; τ_b]
    Af = [zeros(3,3) I zeros(3,3); zeros(3,6) (rad2arcsec/Nm2μNm).*inv(J); zeros(3,9)]
    Bf = Array([zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J); zeros(3,3)]) #control input jacobian
    Cf = [Array(Diagonal(ones(6))) zeros(6,3)]

    #Convert to discrete time
    Hf = exp(dt.*[Af Bf; zeros(3,12)])
    Af = Hf[1:9,1:9]
    Bf = Hf[1:9,10:12]

    #LQR Controller Design with Bias Torque
    #Note that dlqr doesn't work since the system is not strictly controllable
    Q = Array(Diagonal([1.0*ones(3); 10.0*ones(3); zeros(3)])) # NOTE: tune this
    R = Array(Diagonal(0.01.*ones(3)))                         # NOTE: tune this

    P = Q
    K = zeros(3,9)
    for k = 1:100
        K = (R+Bf'*P*Bf)\Bf'*P*Af
        P = Q + K'*R*K + (Af-Bf*K)'*P*(Af-Bf*K)
    end

    #Closed-loop sim
    x = zeros(6,length(t))
    u = zeros(3,length(t)-1)
    # initial condition
    x[:,1] .= [0*randn(3); 0*randn(3)]

    x̄ = zeros(9,length(t)) #Filter state includes torque bias
    x̄[:,1] = [zeros(6);τ[:,1]]
    P = zeros(9,9,length(t))
    P[:,:,1] .= Array(Diagonal([100*ones(6); ones(3)]))

    for k = 1:(length(t)-1)
        #Run controller one step
        u[:,k] .= -K*x̄[:,k] #[x[:,k]; τ[:,k]]

        # true dynamics update
        x[:,k+1] .= A*x[:,k] + B*u[:,k] + G*τ[:,k]

        # true measurement update
        y = C*x[:,k+1] + sqrt(W)*randn(6) #Generate measurement with appropriate noise

        # KF predict
        xp = Af*x̄[:,k] + Bf*u[:,k] #State prediction
        Pp = Af*P[:,:,k]*Af' + V #Prediction covariance

        # KF innovate
        z = y - Cf*xp #Innovation
        S = Cf*Pp*Cf' + W #Innovation covariance
        L = Pp*Cf'*inv(S) #Kalman gain

        # KF update
        x̄[:,k+1] .= xp + L*z #Measurement update

        # covariance update (Joseph Form for psd stabilization)
        P[:,:,k+1] .= (I - L*Cf)*Pp*(I - L*Cf)' + L*W*L'
    end


    # compute RMS pointing error in arcseconds
    Nt = length(t)
    Xvar = 0.0
    for k = 1:Nt
        Xvar += (1/Nt)*(x[1:3,k]'*x[1:3,k])
    end
    RMSerr = sqrt(Xvar)
    return RMSerr

end

function run_mc()
    """Runs a Monte-Carlo of main_runs×46 simulations varying sensing error"""

    # number of runs through the sensing error range
    main_runs = 1000
    rms_mat = zeros(main_runs,46)
    σ_telescope_range = zeros(46)

    # creating an evenly spaced sensing range on a log scale
    σ_telescope_range[1] = 0.0001 #arcsec^2 1-sigma at 1 Hz
    for i = 2:length(σ_telescope_range)
        σ_telescope_range[i] = 1.3*σ_telescope_range[i-1]
    end

    # main loop
    @showprogress "simulating..." for ii = 1:main_runs
        rms_vec = zeros(length(σ_telescope_range))
        for i = 1:length(rms_vec)
            W_telescope_var = σ_telescope_range[i]^2 # var = σ²
            rms_vec[i] = run_single_sim(τ_hist,B_hist_b,J,W_telescope_var)
        end
        rms_mat[ii,:] = rms_vec
    end

    # get mean pointing error and σ bounds
    lower_bound = zeros(46)
    mean_line = zeros(46)
    upper_bound = zeros(46)
    for i = 1:length(mean_line)
        mean_line[i] = mean(rms_mat[:,i])
        σ = std(rms_mat[:,i])
        lower_bound[i] = mean_line[i] - 2*σ
        upper_bound[i] = mean_line[i] + 2*σ
    end

    mat"
    figure
    hold on
    plot($σ_telescope_range,$mean_line,'b')
    plot($σ_telescope_range,$lower_bound,'r')
    plot($σ_telescope_range,$upper_bound,'r')
    xlabel('Telescope Sensing Error (1-sigma, arcseconds)')
    ylabel('RMS Body Pointing Error (arcseconds)')
    xlim([1e-4,10])
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold off
    "
@save "mc_pointing_rms.jld2" σ_telescope_range mean_line lower_bound upper_bound
end

run_mc()
