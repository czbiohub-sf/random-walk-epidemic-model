using Distributions
using Random
using XLSX

function moore_neighbors(site)
    [
        [site[1] + 1, site[2]],
        [site[1] + 1, site[2] + 1],
        [site[1] - 1, site[2] + 1],
        [site[1] - 1, site[2]-1],
        [site[1] - 1, site[2]],
        [site[1], site[2] + 1],
        [site[1], site[2] - 1],
        [site[1] + 1, site[2] - 1]
    ]
end

function probinfected(x, p1, p2)
    if x <= 0
        return p1
    else
        return p2
    end
end

function fill_lattice_data(sheet, arr, savetimes)
    #=
    We have four dimensions in the removed array:
    data[1][2][3][4]
    1: the iteration in n_iter
    2: the timestep in savetimes
    3: the site in the data (not in any particular order except for how it was added chronologically)
    4: the element in the data, e.g. for an infected site there are 4 pieces of data
    [the x coor, y coor, lifetime, and # of ppl inf or something]
    for a removed site data, its just two pieces of data, [x coor, y coor]
    =#

    dim1 = size(arr); dim1 = dim1[1] # This number should be constant
    dim2 = size(arr[1]); dim2 = dim2[1] # This number should be constant
    # dim 3 is not constant - the number of sites in each realization or timestep is not fixed
    # dim 4 is not necessary for our purposes - we just join the array elements into a single string

    for i in 1:dim1
        for j in 1:dim2
            # calculate the value of dim3 since it changes each time
            dim3 = size(arr[i][j]); dim3 = dim3[1]
            XLSX.setdata!(sheet, XLSX.CellRef( 1, (dim2+1)*(i-1) + j ), savetimes[j])
            for k in 1:dim3
                sitedata = arr[i][j][k]
                sitedata = join(sitedata, ",")
                XLSX.setdata!(sheet, XLSX.CellRef( 1 + k, (dim2+1)*(i-1) + j ), sitedata)
            end
        end
    end
end

function simpler_jump(x, y, a = 1)
    theta = rand(Uniform(0, 2 * π))
    u = rand(Uniform(0, 1))

    r = a/(3 * u)^(1/3)

    xfl = r * sin(theta)
    yfl = r * cos(theta)

    (floor(Int, x + xfl + 1/2),
    floor(Int, y + yfl + 1/2), r)
end

function get_cluster_points(start_point, removed)
    neighboringremoved = true
    cluster = Set([start_point])
    cluster_added_old = cluster
    cluster_added_new = Set()
    neighbors_visited = Set()

    #cluster_added_old contains the new sites added to the central cluster in the previous time step
    #cluster_added_new contains the new sites which are currently being added in a given time step
    #each loop, we check the neighbors of cluster_added_old
    #if cluster_added_new is empty, we end the while loop

    while neighboringremoved
        for site in cluster_added_old
            neighbor_sites = moore_neighbors(site)

            for neighbor_site in neighbor_sites
                if (neighbor_site in removed) && !(neighbor_site in neighbors_visited)
                    push!(cluster, neighbor_site)
                    push!(cluster_added_new, neighbor_site)
                    push!(neighbors_visited, neighbor_site)
                end
            end
        end

        if isempty(cluster_added_new) # we have reached limit of central cluster
            neighboringremoved = false
        end

        cluster_added_old = cluster_added_new
        cluster_added_new = Set()
    end

    cluster
end

function surrounded_by_removed(site, removed)
    [site[1]+1, site[2]] in removed &&
        [site[1]+1, site[2]+1] in removed &&
        [site[1], site[2]+1] in removed &&
        [site[1]-1, site[2]+1] in removed &&
        [site[1]-1, site[2]] in removed &&
        [site[1]-1, site[2]-1] in removed &&
        [site[1], site[2]-1] in removed &&
        [site[1]+1, site[2]-1] in removed
end
