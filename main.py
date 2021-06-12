import csv
from time import sleep
import clipboard
import keyword
import pyautogui
import requests

seqs_information_dict = {}

with open('completed_seqs.csv', newline='', encoding="utf-8") as file:
    done_data = csv.reader(file)
    done_seqs = [seq[0].upper() for seq in done_data if seq[0] != 'SEQUENCE']

with open('ignore_these_seqs.txt', newline='', encoding="utf-8") as file:
    ignore = csv.reader(file)
    ignored_list = [seq[0].upper() for seq in ignore]


with open('long_seqs.txt', encoding='utf-8', newline='') as file2:
    reader = csv.reader(file2)
    all_long_seqs = [seq[0] for seq in reader]

with open('all-seq.txt', newline='', encoding="utf-8") as file:
    data = csv.reader(file)
    required_seq = [seq[0].upper() for seq in data if
                    seq[0].upper() not in ignored_list and seq[0].upper() not in done_seqs]

# ------------------------------------------------


def check_connection():
    while True:
        try:
            print(requests.get('https://www.google.com/').status_code)
            break
        except:
            print('no internet')
            sleep(3)
            check_connection()



def is_availabe(move_to, equal_to, sleeep, num_of_repeat):
    # check_connection()
    sleep(sleeep)
    for i in range(num_of_repeat):
        pyautogui.moveTo(*move_to, duration=A)
        pyautogui.tripleClick()
        pyautogui.hotkey('ctrl', 'c')
        right_copy = clipboard.paste()
        if str(equal_to) in right_copy:
            break
        sleep(sleeep)

    return clipboard.paste()


calc_button = (390, 44)
geometry_button = (470, 366)
topolgy_analysis_button = (653, 358)
other_button = (473, 396)
h_bond_button = (661, 402)

def add_new_seq():
    middle_of_screen = (922, 509)
    ungroup = (942, 710)
    # --------------------------------------------------
    sleep(4)
    seqs_information_dict.update({'SEQUENCE': seq})

    clipboard.copy(seq)
    pyautogui.moveTo(*middle_of_screen, duration=A)
    pyautogui.hotkey('ctrl', 'v')

    sleep(B)
    pyautogui.hotkey('ctrl', 'a')
    sleep(B)

    pyautogui.click(button='right')
    pyautogui.moveTo(*ungroup, duration=A)
    pyautogui.click(button='left')
    pyautogui.hotkey('ctrl', 'a')
    sleep(B)
    pyautogui.hotkey('ctrl', '2')


def calc_atom_count():
    elemntal_analysis = (400, 65)
    elemntal_analysis_ok = (848, 719)
    atom_count_copy = (837, 313)
    exit_atom_count_window = (1176, 281)
    exit_atom_count_window2 = (1098, 352)
    # ---------------------------------------

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*elemntal_analysis, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*elemntal_analysis_ok, duration=A)
    pyautogui.click()
    sleep(D)

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*atom_count_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    try:
        seqs_information_dict.update({'atom count': clipboard.paste().split('>')[4].split('<')[0]})
    except:
        seqs_information_dict.update({'atom count': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*exit_atom_count_window, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*exit_atom_count_window2, duration=A)
    pyautogui.click()


def calc_topolgy_analysis():
    topolgy_analysis_ok = (845, 882)
    asymmetric_atom_count_copy = (607, 317)
    rotatable_bond_count_copy = (603, 360)
    ring_count_copy = (701, 464)
    aromatic_ring_count_copy = (606, 405)
    hetero_ring_count_copy = (590, 445)
    exit_topology_analysis = (1319, 262)
    exit_topology_analysis2 = (1198, 182)
    # ----------------------------------------

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*geometry_button, duration=A)
    pyautogui.moveTo(*topolgy_analysis_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*topolgy_analysis_ok, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check = is_availabe(asymmetric_atom_count_copy, 'count =', C, 100)
    try:
        seqs_information_dict.update({'asymmetric atom count': check.split('=')[1]})
    except:
        seqs_information_dict.update({'asymmetric atom count': check})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*rotatable_bond_count_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'rotatable bond count': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*ring_count_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'ring count': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*aromatic_ring_count_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'aromatic ring count': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*hetero_ring_count_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'hetero ring count': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*exit_topology_analysis, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*exit_topology_analysis2, duration=A)
    pyautogui.click()


def calc_geometrical_decription():
    geometrical_decription_button = (662, 402)
    geometrical_decription_ok = (838, 863)
    van_der_Waals_volume_copy = (685, 384)
    minimal_projection_surface_area_copy = (703, 257)
    maximal_projection_surface_area_copy = (705, 281)
    minimal_projection_radius_copy = (721, 305)
    maximal_projection_radius_copy = (724, 325)
    exit_geometrical_decription = (1446, 182)
    exit_geometrical_decription2 = (1198, 203)
    # -----------------------------------------------------
    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*geometry_button, duration=A)
    pyautogui.moveTo(*geometrical_decription_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*geometrical_decription_ok, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check = is_availabe(van_der_Waals_volume_copy, 'van der Waals volume =', C, 100)
    try:
        seqs_information_dict.update({'the van der waals volume': check.split('=')[1]})
    except:
        seqs_information_dict.update({'van der waals surface area': check})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*minimal_projection_surface_area_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'minimal projection surface area': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*maximal_projection_surface_area_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'maximal projection surface area': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*minimal_projection_radius_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'minimal projection radius': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*maximal_projection_radius_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'maximal projection radius': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*exit_geometrical_decription, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*exit_geometrical_decription2, duration=A)
    pyautogui.click()


def calc_molecular_surface_area():
    molecular_surface_area_button = (661, 463)
    van_der_waals_surface_area = (847, 482)
    ok_button = (842, 686)
    Van_der_Waals_surface_area3d_copy = (1028, 313)
    Van_der_Waals_surface_area3d_exit = (1158, 278)
    solvant_accessable_button = (847, 514)
    ASA_copy = (798, 268)
    ASA_exit = (1177, 233)
    full_exit = (1096, 386)
    # ------------------------------------------

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*geometry_button, duration=A)
    pyautogui.moveTo(*topolgy_analysis_button, duration=A)
    pyautogui.moveTo(*molecular_surface_area_button, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*van_der_waals_surface_area, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*ok_button, duration=A)
    pyautogui.click()

    check = is_availabe(Van_der_Waals_surface_area3d_copy, 'Van der Waals surface area (3D) =', C, 100)
    try:
        seqs_information_dict.update({'van der waals surface area': check.split('=')[1]})
    except:
        seqs_information_dict.update({'van der waals surface area': check})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # sleep(E)
    # pyautogui.moveTo(Van_der_Waals_surface_area3d_copy, duration=A)
    # pyautogui.doubleClick()
    # pyautogui.hotkey('ctrl', 'c')
    # seqs_information_dict.update({'van der waals surface area': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*Van_der_Waals_surface_area3d_exit, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*solvant_accessable_button, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*ok_button, duration=A)
    pyautogui.click()
    pyautogui.hotkey('ctrl', 'a')

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check2 = is_availabe(ASA_copy, 'ASA =', C, 100)
    try:
        seqs_information_dict.update({'solvent accessible surface area': check2.split('=')[1]})
    except:
        seqs_information_dict.update({'solvent accessible surface area': check2})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*ASA_exit, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*full_exit, duration=A)
    pyautogui.click()


def calc_polar_surface_are():
    polar_surfave_area_button = (664, 432)
    polar_surfave_area_button_ok = (843, 639)
    polar_surface_area_copy = (928, 313)
    polar_surface_area_exit = (1178, 285)
    full_exit = (1096, 432)

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*geometry_button, duration=A)
    pyautogui.moveTo(*topolgy_analysis_button, duration=A)
    pyautogui.moveTo(*polar_surfave_area_button, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*polar_surfave_area_button_ok, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check = is_availabe(polar_surface_area_copy, 'Polar surface area =', B, 1000)
    try:
        seqs_information_dict.update({'polar surface area': check.split('=')[1]})
    except:
        seqs_information_dict.update({'polar surface area': check})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*polar_surface_area_exit, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*full_exit, duration=A)
    pyautogui.click()


def calc_polarizability():
    charge_button = (475, 233)
    charge_button2 = (641, 240)
    polarizability_button = (641, 276)
    polarizability_ok = (838, 677)
    polarizability_copy = (842, 315)
    polarizability_exit = (1175, 279)
    polarizability_exit2 = (1101, 391)

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*charge_button, duration=A)
    pyautogui.moveTo(*charge_button2, duration=A)
    pyautogui.moveTo(*polarizability_button, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*polarizability_ok, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check = is_availabe(polarizability_copy, 'molecular =', B, 100)
    try:
        seqs_information_dict.update({'polarizability': check.split('=')[1]})
    except:
        seqs_information_dict.update({'polarizability': check})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*polarizability_exit, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*polarizability_exit2, duration=A)
    pyautogui.click()


def calc_H_bond():
    h_bond_ok = (840, 763)
    doner_copy = (1035, 627)
    acceptor_copy = (1130, 669)
    close_window = (1239, 564)
    close_window2 = (1298, 729)
    full_close = (1096, 313)

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*other_button, duration=A)
    pyautogui.moveTo(*h_bond_button, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*h_bond_ok, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check = is_availabe(doner_copy, 'Donor count =', C, 100)
    try:
        seqs_information_dict.update({'H-bond donor count': check.split('=')[1]})
    except:
        seqs_information_dict.update({'H-bond donor count': check})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pyautogui.moveTo(*acceptor_copy, duration=A)
    pyautogui.doubleClick()
    pyautogui.hotkey('ctrl', 'c')
    seqs_information_dict.update({'H-bond acceptor_copy count': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*close_window, duration=A)
    pyautogui.click(button='right')
    pyautogui.moveTo(*close_window2, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*full_close, duration=A)
    pyautogui.click()


def calc_partitioning(seq):
    partitioning_button = (439, 140)
    logP = (608, 138)
    logP_ok = (851, 718)
    logP_copy = (978, 305)
    logP_exit = (1177, 274)
    logP_exit2 = (1113, 352)

    lodD = (611, 174)
    lodD_ok = (842, 763)
    logD_copy = (1008, 708)
    logD_exit = (1263, 568)
    logD_exit2 = (1315, 730)
    logD_exit3 = (1113, 311)

    HLB = (605, 205)
    HLB_ok = (840, 631)
    HLB_copy = (883, 313)
    HLB_exit = (1175, 280)
    HLB_exit2 = (1099, 438)

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*partitioning_button, duration=A)
    pyautogui.moveTo(*logP, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*logP_ok, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    val ='logP of nonionic species ='

    if len(seq)<=4:
        val = 'logP ='

    check = is_availabe(logP_copy, val, C, 100)
    try:
        seqs_information_dict.update({'partition coefficient (logP)': check.split('=')[1]})
    except:
        seqs_information_dict.update({'partition coefficient (logP)': check})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*logP_exit, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*logP_exit2, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*partitioning_button, duration=A)
    pyautogui.moveTo(*logP, duration=A)
    pyautogui.moveTo(*lodD, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*lodD_ok, duration=A)
    pyautogui.click()

    pyautogui.hotkey('ctrl', 'a')
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check2 = is_availabe(logD_copy, '7.40', C, 100)
    try:
        seqs_information_dict.update({'logD': check2.split()[1]})
    except:
        seqs_information_dict.update({'logD': check2})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*logD_exit, duration=A)
    pyautogui.click(button='right')
    pyautogui.moveTo(*logD_exit2, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*logD_exit3, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*partitioning_button, duration=A)
    pyautogui.moveTo(*logP, duration=A)
    pyautogui.moveTo(*HLB, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*HLB_ok, duration=A)
    pyautogui.click()
    pyautogui.hotkey('ctrl', 'a')

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    check3 = is_availabe(HLB_copy, 'Chemaxon HLB =', B, 100)
    try:
        seqs_information_dict.update({'HLB': check3.split('=')[1]})
    except:
        seqs_information_dict.update({'HLB': check3})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*HLB_exit, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*HLB_exit2, duration=A)
    pyautogui.click()


def calc_intrinsic():
    intrinsic = (441, 168)
    intrinsic_button = (656, 170)
    intrinsic_copy = (302, 547)
    full_exit = (1240, 15)

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*intrinsic, duration=A)
    pyautogui.moveTo(*intrinsic_button, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    why = is_availabe(intrinsic_copy, 'Intrinsic solubility:', E, 600)
    try:
        seqs_information_dict.update({'intrinsic solubility': why.split(':')[1].replace('logS', '')})
    except:
        seqs_information_dict.update({'intrinsic solubility': why})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*full_exit, duration=A)
    pyautogui.click()


def calc_refractivity():
    refractivity = (643, 464)
    refractivity_ok = (839, 647)
    refractivity_copy = (877, 313)
    full_exit0 = (1174, 277)
    full_exit1 = (1097, 428)

    pyautogui.hotkey('ctrl', 'a')

    pyautogui.moveTo(*calc_button, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*other_button, duration=A)
    pyautogui.moveTo(*h_bond_button, duration=A)
    pyautogui.moveTo(*refractivity, duration=A)
    pyautogui.click()

    pyautogui.moveTo(*refractivity_ok, duration=A)
    pyautogui.click()

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    is_availabe(refractivity_copy, 'Refractivity =', B, 100)
    try:
        seqs_information_dict.update({'refractivity': clipboard.paste().split('=')[1]})
    except:
        seqs_information_dict.update({'refractivity': clipboard.paste()})
    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    pyautogui.moveTo(*full_exit0, duration=A)
    pyautogui.click()
    pyautogui.moveTo(*full_exit1, duration=A)
    pyautogui.click()

    pyautogui.hotkey('ctrl', 'a')
    pyautogui.press('del')


pyautogui.FAILSAFE = True
# print(keyword.iskeyword('q'))

while keyword.iskeyword('q') != True:
    for seq in required_seq:

        A = 0.1
        B = 0.3
        C = 1
        D = 2
        E = 3
        seq = seq.upper()

        if "X" in seq or '-' in seq:
            continue

        if len(seq) >= 55:
            if seq not in all_long_seqs:
                with open('long_seqs.txt', encoding='utf-8', newline='', mode='a') as file2:
                    write1 = csv.writer(file2)
                    write1.writerow(seq.split())
            continue
            # C += 12
            # D += 12
            # E += 44


        elif len(seq) >= 45:
            A += 0.2
            B += 1
            C += 6
            D += 6
            E += 6

        elif len(seq) >= 30:
            C += 2
            D += 2
            E += 2

        print('Start with this one:', seq)
        print(len(seq))

        add_new_seq()
        calc_atom_count()
        calc_topolgy_analysis()
        calc_geometrical_decription()
        calc_molecular_surface_area()
        calc_polar_surface_are()
        calc_polarizability()
        calc_H_bond()
        calc_partitioning(seq)
        calc_intrinsic()
        calc_refractivity()

        print(seqs_information_dict)

        keys = seqs_information_dict.keys()
        with open('completed_seqs.csv', 'a', newline='', encoding="utf-8")  as file:
            dict_writer = csv.DictWriter(file, keys)
            if file.tell() == 0:
                dict_writer.writeheader()
            dict_writer.writerow(seqs_information_dict)

        print('#' * 10)
        print('Done!:', seq)
        print(len(seq))
        print('#' * 50)

#
# while True:
#     sleep(1)
#     print(pyautogui.position())
#
