#Name: Labeling generic residue numbers assigned via GPCRdb Tools...
#Command: pythonrun maestro_label_generic_numbers.panel

__doc__ = """
$Revision: 0.1 $
$Date: 2015/02/04 $
$Author: Stefan Mordalski $
"""


from PyQt4 import QtCore, QtGui
import schrodinger.ui.qt.appframework as appframework
from schrodinger.maestro import maestro
from schrodinger.structutils import analyze

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class generic_labeling_ui(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(197, 303)
        Dialog.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.gridLayout = QtGui.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.gpcrdb_button = QtGui.QPushButton(Dialog)
        self.gpcrdb_button.setObjectName(_fromUtf8("gpcrdb_button"))
        self.verticalLayout.addWidget(self.gpcrdb_button)
        self.bw_button = QtGui.QPushButton(Dialog)
        self.bw_button.setObjectName(_fromUtf8("bw_button"))
        self.verticalLayout.addWidget(self.bw_button)
        self.both_button = QtGui.QPushButton(Dialog)
        self.both_button.setObjectName(_fromUtf8("both_button"))
        self.verticalLayout.addWidget(self.both_button)
        self.clear_label_button = QtGui.QPushButton(Dialog)
        self.clear_label_button.setObjectName(_fromUtf8("clear_label_button"))
        self.verticalLayout.addWidget(self.clear_label_button)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        

    def retranslateUi(self, Dialog):
        self.gpcrdb_button.setText(_translate("Dialog", "Label GPCRdb", None))
        self.bw_button.setText(_translate("Dialog", "Label B-W", None))
        self.both_button.setText(_translate("Dialog", "Label both", None))
        self.clear_label_button.setText(_translate("Dialog", "Clear labels", None))


class generic_labeling_gui(appframework.AppFramework):
    
    def __init__(self):
        appframework.AppFramework.__init__(
            self,
            ui = generic_labeling_ui(),
            title = "Label generic numbers",
            subwindow = True
        )
        self.ui.gpcrdb_button.clicked.connect(self.label_gpcrdb)
        self.ui.bw_button.clicked.connect(self.label_bw)
        self.ui.both_button.clicked.connect(self.label_both)
        self.ui.clear_label_button.clicked.connect(self.clear_labels)




    def evaluate_workspace(self):
        st = maestro.workspace_get()
        if st.atom_total == 0:
            msg = QtGui.QMessageBox.warning(self, "Empty workspace", "Please place some structure in the workspace!", QtGui.QMessageBox.Ok)
            return 0
        return 1


    def label_gpcrdb(self):
        
        if self.evaluate_workspace():

            maestro.command("labelatom anum=false element=false pdbbfactor=true a. CA and a.pdbbfactor > 0 and a.pdbbfactor < 8.1")

            st = maestro.workspace_get()
            maestro.workspace_set(st, regenerate_markers = True)
            atom_list = analyze.get_atoms_from_asl(st, "a. CA and a.pdbbfactor < 0 and a.pdbbfactor > -8.1")
            for atom in atom_list:
                maestro.command('labelusertextatomupdate utext=%s a. %i' %(((-float(atom.temperature_factor))+0.001), atom.index))


    def label_bw(self):

        if self.evaluate_workspace():
            label_asl = "((atom.ptype N) AND (atom.pdbbfactor >0) AND (atom.pdbbfactor <8.1))"
            maestro.command("labelatom anum=false element=false pdbbfactor=true %s" % label_asl )


    def label_both(self):

        if self.evaluate_workspace():
            
            maestro.command( "labelclear all" ) 

            st = maestro.workspace_get()
            gpcrdb_list = analyze.get_atoms_from_asl(st, "a. CA and a.pdbbfactor > -8.1 and a.pdbbfactor < 8.1")
            bw_list = analyze.get_atoms_from_asl(st, "a.pt N AND a.pdbbfactor > 0 AND a.pdbbfactor < 8.1")
            pairs = zip([x for x in gpcrdb_list], [y for y in bw_list])

            normal = []
            for gpcrdb, bw in pairs:
                if gpcrdb.resnum == bw.resnum and gpcrdb.temperature_factor == bw.temperature_factor:
                    normal.append(str(gpcrdb.index))
                if gpcrdb.resnum == bw.resnum and gpcrdb.temperature_factor != bw.temperature_factor:
                    if gpcrdb.temperature_factor < 0:
                        maestro.command('labelusertextatomupdate utext=%.2fx%s a. %i' %(bw.temperature_factor, ((-float(gpcrdb.temperature_factor))+0.001), gpcrdb.index))
                    else:
                        maestro.command('labelusertextatomupdate utext=%.2fx%.2f a. %i' %(bw.temperature_factor, gpcrdb.temperature_factor, gpcrdb.index))

            maestro.command("labelatom anum=false element=false pdbbfactor=true a. %s" % ','.join(normal) )

    def clear_labels(self):

        maestro.command( "labelclear all" ) 


test_app = None
        
def panel():

    global test_app

    if not test_app:
        test_app = generic_labeling_gui()
    
    test_app.show()


if __name__ == "__main__":

    print "You must run this app from within the Maestro"
