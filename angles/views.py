from django.shortcuts import render
from django.conf import settings
from django.views.generic import TemplateView, View
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.decorators.cache import cache_page

import contactnetwork.pdb as pdb
from structure.models import Structure
from residue.models import Residue
from angles.models import ResidueAngle as Angle

import Bio.PDB
import copy
import io
from collections import OrderedDict
import numpy as np
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d
import freesasa
import scipy.stats as stats

def angleAnalysis(request):
    """
    Show angle analysis site
    """
    return render(request, 'angles/angleanalysis.html')

# TODO: rename and move to tool
class buildAngles(TemplateView):

    template_name = "test.html"
    def get_context_data(self, **kwargs):
        dblist = []
        context = super(testTemplate, self).get_context_data(**kwargs)
        extra_pca = True

        ###########################################################################
        ############################ Helper  Functions ############################
        ###########################################################################

        failed = []

        # Get all structures
        references = Structure.objects.filter(protein_conformation__protein__family__slug__startswith="001").exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')

        # DEBUG for a specific PDB
        #references = Structure.objects.filter(pdb_code__index="2R4R").filter(protein_conformation__protein__family__slug__startswith="001").exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        references = list(references)

        pids = [ref.protein_conformation.protein.id for ref in references]

        qset = Residue.objects.filter(protein_conformation__protein__id__in=pids)
        qset = qset.filter(generic_number__label__regex=r'^[1-7]x[0-9]+').order_by('-protein_conformation__protein','-generic_number__label')
        qset = list(qset.prefetch_related('generic_number', 'protein_conformation__protein','protein_conformation__state'))

        res_dict = {ref.pdb_code.index:qgen(ref.protein_conformation.protein,qset) for ref in references}

        #######################################################################
        ######################### Start of main loop ##########################
        #######################################################################

        angle_dict = [{},{},{},{}]
        median_dict = [{},{},{},{}]

        for reference in references:
            preferred_chain = reference.preferred_chain.split(',')[0]
            pdb_code = reference.pdb_code.index
#            print(pdb_code)

            try:
            #if True:
                structure = load_pdb_var(pdb_code,reference.pdb_data.pdb)
                pchain = structure[0][preferred_chain]
                state_id = reference.protein_conformation.state.id

                #######################################################################
                ###################### prepare and evaluate query #####################

                db_reslist = res_dict[pdb_code]
                #print(db_reslist)

                #######################################################################
                ######################### filter data from db #########################

                def reslist_gen(x):
                    try:
                        while db_reslist[-1].generic_number.label[0] == x:
                            yield db_reslist.pop()
                    except IndexError:
                        pass

                # when gdict is not needed the helper can be removed
                #db_tmlist = [[(' ',r.sequence_number,' ') for r in reslist_gen(x) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
                # db_helper = [[(r,r.sequence_number) for r in reslist_gen(x) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
                db_helper = [[(r,r.sequence_number) for r in reslist_gen(x) if r.sequence_number in pchain] for x in ["1","2","3","4","5","6","7"]]
                gdict = {r[1]:r[0] for hlist in db_helper for r in hlist}
                db_tmlist = [[(' ',r[1],' ') for r in sl] for sl in db_helper]
                db_set = set(db_tmlist[0]+db_tmlist[1]+db_tmlist[2]+db_tmlist[3]+db_tmlist[4]+db_tmlist[5]+db_tmlist[6])


                #######################################################################
                ##################### Angles/dihedrals residues #######################

                polychain = [ residue for residue in pchain if Bio.PDB.Polypeptide.is_aa(residue) and "CA" in residue]
                poly = Bio.PDB.Polypeptide.Polypeptide(polychain)
                poly.get_phi_psi_list() # backbone dihedrals
                poly.get_theta_list() # angle three consecutive Ca atoms
                poly.get_tau_list() # dihedral four consecutive Ca atoms

                # TODO: extend with Chi1-5?
                # https://gist.github.com/lennax/0f5f65ddbfa278713f58
                # Definition http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
                # http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-pysrc.html#Polypeptide.get_phi_psi_list

                dihedrals = {}
                for r in poly:
                  angle_list = ["PHI", "PSI", "THETA", "TAU"]
                  for angle in angle_list:
                      if angle not in r.xtra:
                          r.xtra[angle] = None
                  dihedrals[r.id[1]] = [r.xtra["PHI"], r.xtra["PSI"], r.xtra["THETA"], r.xtra["TAU"]]

                ### clean the structure to solely the 7TM bundle
                recurse(structure, [[0], preferred_chain, db_set])

                # Extra: remove hydrogens from structure (e.g. 5VRA)
                for residue in structure[0][preferred_chain]:
                    for id in [atom.id for atom in residue if atom.element == "H"]:
                        residue.detach_child(id)

                ### AXES through each of the TMs and the TM bundle (center axis)
                hres_list = [np.asarray([pchain[r]["CA"].get_coord() for r in sl], dtype=float) for sl in db_tmlist]
                h_cb_list = [np.asarray([pchain[r]["CB"].get_coord() if "CB" in pchain[r] else cal_pseudo_CB(pchain[r]) for r in sl], dtype=float) for sl in db_tmlist]
                #print(hres_list)
                # fast and fancy way to take the average of N consecutive elements
                N = 3
                hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])
                #print(hres_three)
                ### PCA - determine axis through center + each transmembrane helix
                helix_pcas = [PCA() for i in range(7)]
                helix_pca_vectors = [pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)]

                # Calculate PCA based on the upper (extracellular) half of the GPCR (more stable)
                if extra_pca:
                    helices_mn = np.asarray([np.mean(h, axis=0) for h in hres_three])

                    # Sensitive to length of array, resulting in array when < 6 coordinates
                    # TODO: investigate and make more robust
                    #pos_list = np.asarray([pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]) for i,h in enumerate(hres_three)])

                    pos_list = []
                    for i,h in enumerate(hres_three):
                        if len(h)>6:
                            pos_list.append(pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]))
                        else:
                            pos_list.append(pca_line(PCA(), h))
                    pos_list = np.asarray(pos_list)

                    pos_list = pos_list - (np.mean(pos_list,axis=1)-helices_mn).reshape(-1,1,3)

                    pca = PCA()
                    # TODO store center_vector + vector to reference residue
                    center_vector = pca_line(pca, np.vstack(pos_list))
                else:
                    pca = PCA()
                    center_vector = pca_line(pca, np.vstack(hres_three))

                ### ANGLES
                # Center axis to helix axis to CA
                a_angle = np.concatenate([axes_calc(h,p,pca) for h,p in zip(hres_list,helix_pcas)]).round(3)

                # Center axis to CA to CB
                b_angle = np.concatenate([ca_cb_calc(ca,cb,pca) for ca,cb in zip(hres_list,h_cb_list)]).round(3)


                ### freeSASA (only for TM bundle)
                # SASA calculations - results per atom
                res, trash = freesasa.calcBioPDB(structure)

                # create results dictionary per residue
                asa_list = {}
                atomlist = list(pchain.get_atoms())
                for i in range(res.nAtoms()):
                    resnum = atomlist[i].get_parent().id[1]
                    if resnum not in asa_list:
                        asa_list[resnum] = 0
                    asa_list[resnum] += res.atomArea(i)

                ### Half-sphere exposure (HSE)
                hse = pdb.HSExposure.HSExposureCB(structure[0])
                hselist = dict([ (x[0].id[1], x[1][1]) if x[1][1] > 0 else 0 for x in hse ])
                continue

                ### Collect all data in database list
                if len(pchain) != len(hselist):
                    raise Exception("\033[91mLength mismatch hse " + pdb_code + "\033[0m")

                if len(pchain) != len(asa_list):
                    raise Exception("\033[91mLength mismatch sasa " + pdb_code + "\033[0m")

                if np.isnan(np.sum(asa_list)):
                    raise Exception("\033[91mNAN sasa " + pdb_code + "\033[0m")

                if np.isnan(np.sum(hselist)):
                    raise Exception("\033[91mNAN hse " + pdb_code + "\033[0m")

                for res,a1,a2,asa,hse in zip(pchain,a_angle,b_angle,asa_list,hselist):
                    dblist.append([gdict[res.id[1]], a1, a2, reference,state_id-1,asa,hse])
                    if gdict[res.id[1]].generic_number.label not in angle_dict[state_id-1]:
                        angle_dict[state_id-1][gdict[res.id[1]].generic_number.label] = [round(a1,3)]
                    else:
                        angle_dict[state_id-1][gdict[res.id[1]].generic_number.label].append(round(a1,3))

            except Exception as e:
            #else:
                print(pdb_code, " - ERROR - ", e)
                failed.append(pdb_code)
                continue

        for i in range(4):
            for key in angle_dict[i]:
                sortlist = np.array(angle_dict[i][key])
                median_dict[i][key] = np.median(sortlist)

        for i, res in enumerate(dblist):
            g = res[0]
            a = res[1]

            templist = copy.copy(angle_dict[res[4]][g.generic_number.label])
            del templist[templist.index(a)]

            std_test = abs(np.average(templist) - int(a))/np.std(templist)
            std_len  = len(templist) - 1
            std = stats.t.cdf(std_test, df=std_len)
            dblist[i].append(0.501 if np.isnan(std) else std)



        dblist = [Angle(residue=g, diff_med=round(abs(median_dict[i][g.generic_number.label]-a1),3), angle=a1, b_angle=a2, structure=ref, sasa=round(asa,3), hse=hse, sign_med=round(sig,3)) for g,a1,a2,ref,i,asa,hse,sig in dblist]
        print("created list")
        print(len(dblist))

        # Store the results
        # faster than updating: deleting and recreating
        Angle.objects.all().delete()
        Angle.objects.bulk_create(dblist,batch_size=5000)

        return context

def get_angles(request):
    data = {error: 0}

    # Request selection
    try:
        pdbs = request.GET.getlist('pdbs[]')
        pdbs = [pdb.upper() for pdb in pdbs]

        # Grab PDB data for (first?) PDB
        query = Angle.objects.filter(structure__pdb_code__index=pdbs[0]).prefetch_related("residue__generic_number")

        # Return prep data
        data['data'] = [[q.residue.generic_number.label,q.residue.sequence_number, q.angle, q.diff_med, q.sign_med, q.hse, q.sasa] for q in query]
    except IndexError:
        data['error'] = 1
        data['errorMessage'] = "No PDB(s) selection provided"

    return JsonResponse(data)

def ServePDB(request, pdbname):
    structure=Structure.objects.filter(pdb_code__index=pdbname.upper())
    if structure.exists():
        structure=structure.get()
    else:
        quit()

    if structure.pdb_data is None:
        quit()

    only_gns = list(structure.protein_conformation.residue_set.exclude(generic_number=None).values_list('protein_segment__slug','sequence_number','generic_number__label').all())
    only_gn = []
    gn_map = []
    segments = {}
    for gn in only_gns:
        only_gn.append(gn[1])
        gn_map.append(gn[2])
        if gn[0] not in segments:
            segments[gn[0]] = []
        segments[gn[0]].append(gn[1])
    data = {}
    data['pdb'] = structure.pdb_data.pdb
    data['only_gn'] = only_gn
    data['gn_map'] = gn_map
    data['segments'] = segments
    data['chain'] = structure.preferred_chain

    return JsonResponse(data)
