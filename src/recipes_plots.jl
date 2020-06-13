
using Plots, UnPack, RollingFunctions

"""
    p1,p2,p3=plot_final(outputs,parameters)
"""
function plot_final(outputs, parameters)

    @unpack U0, UA, THETA, SST, TO, ustar, psimh, psixh = outputs
    @unpack xi, rd, rh, re, SW, LW, SH, LH = outputs
    @unpack dt, n, ndays, cp, dz, rhow0 = parameters

    k1 = plot(
        rollmean(U0[1, :], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        ylabel = "U0",
    )
    k2 = plot(
        rollmean(UA[:, end], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        ylabel = "Utop",
    )
    k3 = plot(
        rollmean(UA[:, end] - U0[1, :], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        ylabel = "Utop-U0",
    )
    k4 = plot(
        rollmean(UA[:, 1], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        ylabel = "U",
    )
    k5 = plot(
        rollmean(THETA[:, 2], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        ylabel = "THETA",
    )
    k6 = plot(
        rollmean(SST, Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        ylabel = "SST",
    )
    k7 = plot(
        rollmean(TO[:, 10], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        ylabel = "T",
    )
    #tt=string("DZ=",dz)
    plt1 = plot(k1, k2, k3, k4, k5, k6, k7, layout = (7, 1), legend = false)

    k1 = plot(
        rollmean(xi[1, :], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "xi",
    )
    k2 = plot(
        rollmean(ustar[1, :], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "ustar",
    )
    k3 = plot(
        rollmean(psimh[1, :], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "psimh",
    )
    k4 = plot(
        rollmean(psixh[1, :], Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "psixh",
    )
    k5 = plot(
        rollmean((rd[1, :] .^ 2), Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "CD",
    )
    k6 = plot(
        rollmean((rh[1, :] .* rd[1, :]), Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "CH",
    )
    k7 = plot(
        rollmean((re[1, :] .* rd[1, :]), Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "CE",
    )
    #tt=string("DZ=",dz)
    plt2 = plot(k1, k2, k3, k4, k5, k6, k7, layout = (7, 1), legend = false)

    k1 = plot(
        rollmean(SW[:, 1] / cp / (dz * rhow0) * dt, Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "SW",
    )
    k2 = plot(
        rollmean(LW[1, :] / cp / (dz * rhow0) * dt, Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "LW",
    )
    k3 = plot(
        rollmean(SH[1, :] / cp / (dz * rhow0) * dt, Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "SH",
    )
    k4 = plot(
        rollmean(LH[1, :] / cp / (dz * rhow0) * dt, Int(24 * 60 * 60 / dt)),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "LH",
    )
    k5 = plot(
        rollmean(
            (SW[:, 1] - LW[1, :] - SH[1, :] - LH[1, :]) / cp / (dz * rhow0) * dt,
            Int(24 * 60 * 60 / dt),
        ),
        xticks = (0:24*60*60/dt:n, string.(0:ndays)),
        title = "QNET",
    )
    plt3 = plot(k1, k2, k3, k4, k5, layout = (5, 1), legend = false)
    return plt1, plt2, plt3
end

function plot_state()
    f1 = plot(-ZO[1:end-1], TO[i, :], ylabel = ("TO"))#,ylims=(19,21)
    f2 = plot(-ZO[1:end-1], UO[i, :], ylabel = ("UO"))
    f3 = plot(-ZO[1:end-1], kom, ylabel = ("KOm"))
    f4 = plot(ZA, UA[i, :], ylabel = ("UA"))
    f5 = plot(ZA, Tv .- 273.15, ylabel = ("Tv")) #ylim([19 21]);...
    f6 = plot(ZA, qA[i, :], ylabel = ("q"))
    f7 = plot(ZA, kam, ylabel = ("KAm"))
    plt =
        plot(f1, f2, f3, f4, f5, f6, f7, layout = (7, 1), legend = false, titlefontsize = 6)
    display(plt)
end
