# outbreak simulation which records the radius of gyration, cluster mass, and
# cluster surface points

# user inputs at bottom of program

using JLD
include("utils.jl")
using FFTW
using LinearAlgebra
using Dates

function simpler_jump(x, y, a = 1)
    theta = rand(Uniform(0, 2 * π))
    u = rand(Uniform(0, 1))

    r = a/(3 * u)^(1/3)

    xfl = r * sin(theta)
    yfl = r * cos(theta)

    (floor(Int, x + xfl + 1/2),
    floor(Int, y + yfl + 1/2), r)
end

function createpaddedgrid(infectors, removed, visited, gridsize)
    xoffset = yoffset = round(Int, gridsize/2)

    grid = zeros(gridsize, gridsize)

    infector_coords = map(x -> x[1:2], infectors)
    removed_coords = map(x -> x[1:2], removed)
    visited_coords = map(x -> x[1:2], visited)

    for (x, y) in visited_coords
        grid[x + xoffset, y + yoffset] = 0
    end

    for (x, y) in removed_coords
        grid[x + xoffset, y + yoffset] = 1
    end

    for (x, y) in infector_coords
        grid[x + xoffset, y + yoffset] = 1
    end

    grid
end

function outbreakwithjumps_varp(maxtime :: Int,
                           maxtau :: Int,
                           nextstepfcn :: Function,
                           probinfected :: Function,
                           p_arr :: Array{Float64, 1},
                           savetimes :: Array{Int, 1},
                           gridsize,

                           initialinfectors = [[0, 0, 0, 0]],
                           addinfectors = true) :: Tuple{
                                                   Array{Array{Array{Int64,N} where N,1},1}, #infectorsTosave
                                                   Array{Array{Array{Int64,N} where N,1},1}, #removedTosave
                                                   Array{Array{Array{Int64,N} where N,1},1}, #visitedTosave
                                                   Array{Float64, 1}, #radiusofgyr_inf
                                                   Array{Float64, 1}, #radiusofgyr_all
                                                   Array{Float64, 1}, #radiusofgyr_cluster
                                                   Array{Int64, 1}, #surfacepts
                                                   Array{Array{Float64,1},1}, #corr_length
                                                   Array{Int64, 1} #centralclustermass
                                                   }
    # Tracking the different agents
    x_start = 0
    y_start = 0

    newinfections = Array{Int64, 1}()
    infectors = initialinfectors
    numinfected = [length(infectors)]
    removed = Set([[x_start, y_start]])
    visited = Set([[x_start, y_start]])
    numvisited = zeros(Int, maxtime)
    numsusceptiblesvisited = zeros(Int, maxtime)
    r_inf = Array{Tuple{Int64, Int64}, 1}()
    graphcomponents = Array{Any, 1}()

    interior = Set()

    radiusofgyr_inf = zeros(Float64, length(savetimes))
    radiusofgyr_all = zeros(Float64, length(savetimes))

    com_inf_x = zeros(Float64, length(savetimes))
    com_inf_y = zeros(Float64, length(savetimes))
    com_all_x = zeros(Float64, length(savetimes))
    com_all_y = zeros(Float64, length(savetimes))

    # since cluster is just removed (not inf), there is only one rofgyr
    radiusofgyr_cluster = zeros(Float64, length(savetimes))

    com_cluster_x = zeros(Float64, length(savetimes))
    com_cluster_y = zeros(Float64, length(savetimes))

    surfacepts = zeros(Int, length(savetimes))
    centralclustermass = zeros(Int, length(savetimes))
    numremoved = zeros(Int, maxtime)

    corr_length = fill(Float64[], length(savetimes))

    removedTosave = fill(Array{Int}[], length(savetimes))
    visitedTosave = fill(Array{Int}[], length(savetimes))
    infectorsTosave = fill(Array{Int}[], length(savetimes))

    savetimeidx = 1

    for time in 1:maxtime
        ninfectors = length(infectors)
        ninfected = 0
        newinfections_now = 0

        if time != 1
            numvisited[time] = numvisited[time - 1]
            numsusceptiblesvisited[time] = numsusceptiblesvisited[time - 1]
        end

        for i in 1:ninfectors
            infector = infectors[i]
            (infector[1], infector[2]) = simpler_jump(infector[1], infector[2])
            infector[3] += 1

            inremoved = infector[1:2] in removed
            alreadyvisited = infector[1:2] in visited

            if !alreadyvisited
                push!(visited, infector[1:2])
                numvisited[time] += 1
            end

            isinfecting = rand() <= probinfected(infector[1], p_arr[1], p_arr[2])

            if !alreadyvisited && !isinfecting
                numsusceptiblesvisited[time] += 1
            end

            if !inremoved && isinfecting && alreadyvisited
                numsusceptiblesvisited[time] -= 1
            end

            if !inremoved && isinfecting
                push!(removed, infector[1:2])
            end

            if !inremoved && isinfecting
                newinfections_now += 1
                if !addinfectors
                    ninfected += 1
                else
                    push!(infectors, [infector[1:2]..., 0, 0])
            end

                infector[4] += 1
            end
        end

        if time in savetimes
            i = findall(x->x==time, savetimes)
            i = i[1]

            if !isempty(infectors)
                com_inf_x[i] = sum(infectors)[1]/length(infectors)
                com_inf_y[i] = sum(infectors)[2]/length(infectors)
                com_all_x[i] = (sum(infectors)[1] + sum(removed)[1]) / (length(infectors) + length(removed))
                com_all_y[i] = (sum(infectors)[2] + sum(removed)[2]) / (length(infectors) + length(removed))
            else
                com_inf_x[i] = NaN
                com_inf_y[i] = NaN
                com_all_x[i] = sum(removed)[1] / length(removed)
                com_all_y[i] = sum(removed)[2] / length(removed)
            end

            for site in infectors
                radiusofgyr_inf[i] += (site[1] - com_inf_x[i])^2 + (site[2] - com_inf_y[i])^2
                radiusofgyr_all[i] += (site[1] - com_all_x[i])^2 + (site[2] - com_all_y[i])^2
            end

            for site in removed
                radiusofgyr_all[i] += (site[1] - com_all_x[i][1])^2 + (site[2] - com_all_y[i][1])^2
            end

            if ninfectors != 0
                radiusofgyr_inf[i] = sqrt(radiusofgyr_inf[i]/ninfectors)
            else
                radiusofgyr_inf[i] = radiusofgyr_inf[i-1]
            end
            radiusofgyr_all[i] = sqrt(radiusofgyr_all[i]/(ninfectors+length(removed)))
        end

        push!(newinfections, newinfections_now)

        infectedtoremove = filter(x -> x[3] >= maxtau, infectors)
        for infector in infectedtoremove
            push!(r_inf, (infector[4], time))
        end

        filter!(x -> x[3] < maxtau, infectors)

        if addinfectors
            push!(numinfected, length(infectors))
        else
            push!(numinfected, ninfected)
        end

        if time in savetimes
            idx = findall(x -> x == time, savetimes)
            idx = idx[1]

            for site in setdiff(removed, interior)
                onsurface = true

                if surrounded_by_removed(site, removed)
                    onsurface = false
                    push!(interior, site)
                end

                if onsurface
                    surfacepts[idx] += 1
                end
            end

            cluster = get_cluster_points([x_start, y_start], removed)
            clustermass = length(cluster)
            centralclustermass[idx] = clustermass

            com_cluster_x[idx] = sum(cluster)[1]/clustermass
            com_cluster_y[idx] = sum(cluster)[2]/clustermass

            for site in cluster
                radiusofgyr_cluster[idx] += (site[1] - com_cluster_x[i])^2 + (site[2] - com_cluster_y[i])^2
            end
            radiusofgyr_cluster[idx] = sqrt(radiusofgyr_cluster[idx]/clustermass)

            # grid data for FFT
            data = createpaddedgrid(infectors, collect(removed), collect(visited), gridsize)
            data = complex(data) # Since output matrix is complex (FFT: R -> C), input needs to be as well

            L = size(data, 1) # data should be LxL matrix

            # pre-allocated complex matrix that is same size of data
            phi = Array{Complex{Float64},2}(undef, L, L)
            C = Array{Float64,2}(undef, L, L)

            # Planning the FFTs
            # See documentation:
            # FFT: https://juliamath.github.io/AbstractFFTs.jl/stable/api/#Public-Interface-1
            # R2R FFT: https://juliamath.github.io/FFTW.jl/latest/fft.html
            pdir = FFTW.plan_fft(data, [1, 2])
            pinv = FFTW.plan_r2r(C, [1, 2])
            # already defined plans as input to the outbreak simulation
            # FFTfwdplan, FFTinvplan

            # FFT
            LinearAlgebra.mul!(phi, pdir, data) # do not normalize this result. the pinv plan will account for this.

            # Square modulus
            F = real(phi).^2 + imag(phi).^2

            # Inverse FFT
            LinearAlgebra.ldiv!(C, pinv, F)

            C_onedim = (C[:, 1] + C[1, :])/2

            corr_length[idx] = C_onedim
        end

        numremoved[time] = length(removed)

        if time in savetimes
            removedTosave[savetimeidx] = collect(removed)
            visitedTosave[savetimeidx] = collect(visited)
            infectorsTosave[savetimeidx] = infectors
            savetimeidx = savetimeidx + 1
        end

        currenttime = Dates.Time(Dates.now()); currentdate = Dates.today()
        #println(string("Finished simulation for iteration ", iternum, " timestep ", time, " at ", currenttime, " on ", currentdate, "."))
    end


    removed = collect(removed)
    visited = collect(visited)

    (infectorsTosave, removedTosave, visitedTosave,
    radiusofgyr_inf, radiusofgyr_all, radiusofgyr_cluster,
    surfacepts, corr_length, centralclustermass)
end

function collectandexportdata(maxtime, p, tau, exportfilename, n_iter, savetimes, gridsize) # , FFTfwdplan, FFTinvplan
    p_arr = [p, p]
    counter = Threads.Atomic{Int}(0)

    listofradiusofgyr_inf = fill(Float64[], n_iter)
    listofradiusofgyr_all = fill(Float64[], n_iter)
    listofradiusofgyr_cluster = fill(Float64[], n_iter)
    listofpercent = fill(Float64[], n_iter)

    # number of removed, inf, and visited
    listofremoved = zeros(Int, n_iter)
    listofinfected = zeros(Int, n_iter)
    listofvisited = zeros(Int, n_iter)

    # to access data in the arrays above, you need 4 indices
    # data[1][2][3][4]
    # 1: the iteration in n_iter
    # 2: the timestep in savetimes
    # 3: the site in the data (not in any particular order except for how it was added chronologically)
    # 4: the element in the data, e.g. for an infected site there are 4 pieces of data
    # [the x coor, y coor, lifetime, and # of ppl inf or something]
    # for a removed site data, its just two pieces of data, [x coor, y coor]

    listofsurfacepts = fill(Float64[], n_iter)
    listofcorrlength = fill(Array{Float64}[], n_iter)
    listofclustermass = fill(Int[], n_iter)

    currenttime = Dates.Time(Dates.now()); currentdate = Dates.today()
    println(string("Beginning simulation for ", exportfilename, " at ", currenttime, " on ", currentdate, "."))

    # counter of number of iterations completed
    # for some reason initializing it at 1 actually sets it to zero?
    x = Threads.Atomic{Int}(1)

    t1 = time()
    t_before = t1
    @Threads.threads for i in 1:n_iter
        # (infections, infectors, removed, graphs) = outbreakwithjumps(time, τ, simpler_jump, times, p)
        currenttime = Dates.Time(Dates.now()); currentdate = Dates.today()
        #println(string("Beginning iteration number ", i, " at ", currenttime, " on ", currentdate, "."))

        (infectorsTosave, removedTosave, visitedTosave,
            radiusofgyr_inf, radiusofgyr_all, radiusofgyr_cluster,
            surfacepts, corr_length, clustermass) =
            outbreakwithjumps_varp(maxtime, tau, simpler_jump, probinfected, p_arr, savetimes, gridsize, i)

        currenttime = Dates.Time(Dates.now()); currentdate = Dates.today()
        #println(string("Finished iteration number ", i, " at ", currenttime, " on ", currentdate, "."))

        listofradiusofgyr_inf[i] = radiusofgyr_inf
        listofradiusofgyr_all[i] = radiusofgyr_all
        listofradiusofgyr_cluster[i] = radiusofgyr_cluster

        listofremoved[i] = length(removedTosave)
        listofinfected[i] = length(infectorsTosave)
        listofvisited[i] = length(visitedTosave)

        listofsurfacepts[i] = surfacepts
        listofcorrlength[i] = corr_length
        listofclustermass[i] = clustermass

        attackrate = Array{Float64, 1}()

        for i = 1:length(savetimes)
            push!(attackrate, ( length(infectorsTosave[i])+length(removedTosave[i]) )/( length(visitedTosave[i]) ) )
        end
        listofpercent[i] = attackrate

        Threads.atomic_add!(x, 1)
        if mod(x[], 10) == 0
            println()
            println(">---------------------------------------------------------------<")
            t2 = time()
            t_elapsed = t2-t1
            currenttime = Dates.Time(Dates.now()); currentdate = Dates.today()
            println(string("Completed ", x[], " iterations in ", round(t_elapsed, digits = 5), " seconds at ", currenttime, " on ", currentdate, "."))
            println(">---------------------------------------------------------------<")
            println()
        end

        #if mod(length(counter), 25) == 0
        #    println(string("Completed ", length(counter), " iterations."))
        #end
    end

    t2 = time()
    t_elapsed = t2-t1
    currenttime = Dates.Time(Dates.now()); currentdate = Dates.today()
    println(string("Simulation for ", exportfilename, " was completed in a total of ", round(t_elapsed, digits = 5), " seconds at ", currenttime, " on ", currentdate, "."))
    println("Beginning to average and save data...")

    percent = mean(listofpercent)
    stddev = std(listofpercent)

    radiusofgyr_inf = mean(listofradiusofgyr_inf)
    radiusofgyr_inf_err = std(listofradiusofgyr_inf)/sqrt(n_iter)
    radiusofgyr_all = mean(listofradiusofgyr_all)
    radiusofgyr_all_err = std(listofradiusofgyr_all)/sqrt(n_iter)
    radiusofgyr_cluster = mean(listofradiusofgyr_cluster)
    radiusofgyr_cluster_err = std(listofradiusofgyr_cluster)/sqrt(n_iter)

    surfacepts = mean(listofsurfacepts)
    surfacepts_err = std(listofsurfacepts)/sqrt(n_iter)
    clustermass = mean(listofclustermass)
    clustermass_err = std(listofclustermass)/sqrt(n_iter)

    numremoved = mean(listofremoved)
    numremoved_err = std(listofremoved)/sqrt(n_iter)
    numinfected = mean(listofinfected)
    numinfected_err = std(listofinfected)/sqrt(n_iter)
    numvisited = mean(listofvisited)
    numvisited_err = std(listofvisited)/sqrt(n_iter)

    fn = string(exportfilename, ".jld")
    save(fn,
        "listofradiusofgyr_inf", listofradiusofgyr_inf, "radiusofgyr_inf", radiusofgyr_inf, "radiusofgyr_inf_err", radiusofgyr_inf_err,
        "listofradiusofgyr_all", listofradiusofgyr_all, "radiusofgyr_all", radiusofgyr_all, "radiusofgyr_all_err", radiusofgyr_all_err,
        "listofradiusofgyr_cluster", listofradiusofgyr_cluster, "radiusofgyr_cluster", radiusofgyr_cluster, "radiusofgyr_cluster_err", radiusofgyr_cluster_err,
        "τ", tau, "p_arr", p_arr, "n_iter", n_iter, "maxtime", maxtime,
        "t_elapsed", t_elapsed,
        "percent", percent, "stddev", stddev,
        "listofpercent", listofpercent,
        "listofsurfacepts", listofsurfacepts, "surfacepts", surfacepts, "surfacepts_err", surfacepts_err,
        "listofclustermass", listofclustermass, "clustermass", clustermass, "clustermass_err", clustermass_err,
        "numremoved", numremoved, "numremoved_err", numremoved_err,
        "numinfected", numinfected, "numinfected_err", numinfected_err,
        "numvisited", numvisited, "numvisited_err", numvisited_err,
        "listofcorrlength", listofcorrlength
        )

    t3 = time()
    currenttime = Dates.Time(Dates.now()); currentdate = Dates.today()
    println(string("Averaged and saved full data as ", fn, " in ", round(t3 - t2, digits = 5), " seconds at ", currenttime, " on ", currentdate, "."))
end

function probinfected(x, p1, p2)
    if x <= 0
        return p1
    else
        return p2
    end
end

# ------------------------------- EDIT INPUTS BELOW ---------------------------

maxtime = 1448
# p = 0.48

n_iter = 500

gridsize = 2750

# rounded powers of sqrt(2) up to 2^(10.5)
savetimes = [1, 2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64, 91, 128, 181, 256, 362, 512, 724, 1024, 1448]

p_inputs = [0.030]
p_names = ["0pt030"]
tau_inputs = [50]
tau_names = ["50"]

for i in 1:length(p_inputs)
    for j in 1:length(tau_inputs)
        p = p_inputs[i]
        p_name = p_names[i]
        tau = tau_inputs[j]
        tau_name = tau_names[j]
        exportfilename = string("fractaldim_tau", tau_name, "_p", p_name)
        collectandexportdata(maxtime, p, tau, exportfilename, n_iter, savetimes, gridsize)
    end
end
