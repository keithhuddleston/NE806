    # if plot:
    #     plt.rcParams.update({'font.size': 14})
    #     fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
    #     ax1.set_xscale("log")
    #     ax1.set_yscale("log")
    #     plt.xlabel('Energy (eV)'), plt.ylabel('Cross Section (Barns)')
    #     ax1.plot(e_vals, xs_vals, color='black')
    #     plt.xlim(min(group_structure), max(group_structure)) 
    #     for k in range(len(group_structure)-1):
    #         ax1.plot([group_structure[k], group_structure[k+1]], 
    #                   [y_vals_group[k], y_vals_group[k]], c='r')
    #         for k in range(len(group_structure)-2):
    #             ax1.plot([group_structure[k+1], group_structure[k+1]], 
    #                       [y_vals_group[k], y_vals_group[k+1]], c='r')
    #     ax2 = ax1.twinx()
    #     ax2.loglog(e_vals, phi)
    #     ax2.set_ylabel('Flux Shape Function Arb. Units')
    #     plt.title(Name+' '+Type+' Cross-Section Evaluated at '
    #               +str(Temp)+r'$\degree\ K\ for\ $'+r'$\sigma_b=$'
    #               +str(xs_dilution)+' Barnes')