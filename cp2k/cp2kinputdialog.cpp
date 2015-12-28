/**********************************************************************
  Cp2kInpuDialog - CP2k Input Dialog

  Copyright (C) Aoyama Iwao

  This file is not yet a part of the Avogadro molecular editor project.
  For more information about Avogadro, see <http://avogadro.cc/>

  Some code is based on Open Babel
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#include "cp2kinputdialog.h"

#include <avogadro/molecule.h>
#include <avogadro/atom.h>

// For use PrimitiveList  GLWidget::selectedPrimitives ()  (Atom Selection)
#include <avogadro/primitive.h>
#include <avogadro/primitivelist.h>
#include <avogadro/glwidget.h>

#include <openbabel/mol.h>

#include <QString>
#include <QTextStream>
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>

using namespace OpenBabel;

namespace Avogadro
{
  Cp2kInputDialog::Cp2kInputDialog(QWidget *parent, Qt::WindowFlags f) : 
                QDialog(parent, f), m_molecule(0), /*m_constraints(0),*/
                m_multiplicity(1), m_charge(0), m_savePath("")
  {
	  ui.setupUi(this);

	  // Connect Qt signal/slots
      // Buttons
	  connect(ui.closeButton, SIGNAL(clicked()), this, SLOT(closeClicked()) );
	  connect(ui.generateButton, SIGNAL(clicked()), this, SLOT(generateClicked()) );
	  connect(ui.resetButton, SIGNAL(clicked()), this, SLOT(resetClicked()) );

	  // Basic tab
	  connect(ui.projectNameLine, SIGNAL(editingFinished()),this, SLOT(setProjectName()) );
	  connect(ui.runTypeCombo, SIGNAL(currentIndexChanged(int)),this, SLOT(setRunType(int)));

	  connect(ui.mmRadioButton, SIGNAL(clicked()), this, SLOT(mmRadioChecked()) );
	  connect(ui.qmRadioButton, SIGNAL(clicked()), this, SLOT(qmRadioChecked()) );
	  connect(ui.qmmmRadioButton, SIGNAL(clicked()), this, SLOT(qmmmRadioChecked()) );

	  // MM tab
	  connect(ui.emaxSplineDoubleSpin, SIGNAL(valueChanged(double)), this, SLOT(setEmaxSpline(double)));
	  connect(ui.ewaldTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(setEwaldType(int)));

	  // QM tab
	  connect(ui.qmMethodCombo,SIGNAL(currentIndexChanged(int)),this, SLOT(setQmMethod(int)) );
	  connect(ui.scfGuessCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(setSCFGuess(int)) );
	  connect(ui.chargeSpin, SIGNAL(valueChanged(int)), this, SLOT(setCharge(int)));
	  connect(ui.multiplicitySpin, SIGNAL(valueChanged(int)), this, SLOT(setMultiplicity(int)));
	  connect(ui.maxSCFSpin, SIGNAL(valueChanged(int)), this, SLOT(setMaxSCF(int)));

	  // DFT tab
	  connect(ui.basisSetCombo,SIGNAL(currentIndexChanged(int)),this, SLOT(setBasisSet(int)) );
	  connect(ui.functionalCombo,SIGNAL(currentIndexChanged(int)),this, SLOT(setFunctional(int)) );
	  connect(ui.nMGridSpin, SIGNAL(valueChanged(int)), this, SLOT(setNMultiGrid(int)));
	  connect(ui.cutOffSpin, SIGNAL(valueChanged(int)), this, SLOT(setCutOff(int)));

	  // DFTB tab
	  connect(ui.sccTypeCombo, SIGNAL(currentIndexChanged(int)),this, SLOT(setDftbSCCType(int)) );
	  connect(ui.dispersionCheck, SIGNAL(stateChanged(int)), this, SLOT(setDftbDisp(int)));

	  // SE tab
	  connect(ui.seMethodCombo, SIGNAL(currentIndexChanged(int)),this, SLOT(setSEMethod(int)) );

	  // QMMM tab
      connect(ui.qmmmElectroCoupleCombo, SIGNAL(currentIndexChanged(int)),this, SLOT(setQMMMElectrostaticCoupling(int)) );

	  // MD tab
	  connect(ui.ensembleCombo, SIGNAL(currentIndexChanged(int)),this, SLOT(setEnsemble(int)) );

	  // Initial settings
	  // Basic tab

	  // MM tab
	  ui.emaxSplineDoubleSpin->setRange( 0.0, 1000.0 );
	  ui.emaxSplineDoubleSpin->setSingleStep( 0.1 );

	  // QM tab
	  ui.chargeSpin->setRange(-10, 10);
	  ui.chargeSpin->setSingleStep(1);

	  ui.multiplicitySpin->setRange(1, 10);
	  ui.multiplicitySpin->setSingleStep(1);

	  ui.maxSCFSpin->setRange( 0, 10000 );
	  ui.maxSCFSpin->setSingleStep(10);

	  // DFT tab
	  ui.nMGridSpin->setRange(1,10);
	  ui.nMGridSpin->setSingleStep(1);

	  ui.cutOffSpin->setRange(1, 1000);
	  ui.cutOffSpin->setSingleStep(1);

	  // DFTB tab

	  // SE tab

	  // QMMM tab

	  // MD tab

	  QSettings settings;
      readSettings(settings);

	  updatePreviewText();
  }

  Cp2kInputDialog::~Cp2kInputDialog()
  {
	  // the following codes apparently does not work...
	  QSettings settings;
      writeSettings(settings);
  }

  void Cp2kInputDialog::setMolecule(Molecule *molecule)
  {
	// Disconnect the old molecule first...
	if (m_molecule) disconnect(m_molecule, 0, this, 0);
   
    m_molecule = molecule;

    // Update the preview text whenever atoms are changed
    connect(m_molecule, SIGNAL(atomRemoved(Atom *)), this, SLOT(updatePreviewText()) );
    connect(m_molecule, SIGNAL(atomAdded(Atom *)), this, SLOT(updatePreviewText()) );
    connect(m_molecule, SIGNAL(atomUpdated(Atom *)), this, SLOT(updatePreviewText()) );
	
    updatePreviewText();

  }


  void Cp2kInputDialog::writeSettings(QSettings &settings) const
  {
    // Basic tab
	 settings.setValue("CP2K/ProjectName", ui.projectNameLine->displayText() );
     settings.setValue("CP2K/RunType", ui.runTypeCombo->currentIndex());
	 settings.setValue("CP2K/MMRadio", m_mmRadioChecked );
	 settings.setValue("CP2K/QMRadio", m_qmRadioChecked );
	 settings.setValue("CP2K/QMMMRadio", m_qmmmRadioChecked );

	 // MM tab
     settings.setValue("CP2K/EmaxSpline", ui.emaxSplineDoubleSpin->value() );
	 settings.setValue("CP2K/EwaldType", ui.ewaldTypeCombo->currentIndex() );

	 // QM tab
	 settings.setValue("CP2K/QmMethod", ui.qmMethodCombo->currentIndex() );
	 settings.setValue("CP2K/SCFGuess", ui.scfGuessCombo->currentIndex() );
	 settings.setValue("CP2K/Charge", ui.chargeSpin->value() );
	 settings.setValue("CP2K/Multiplicity", ui.multiplicitySpin->value() );
	 settings.setValue("CP2K/MaxSCF", ui.maxSCFSpin->value() );

	 // DFT tab
	 settings.setValue( "CP2K/BasisSet", ui.basisSetCombo->currentIndex() );
	 settings.setValue( "CP2K/Functional", ui.functionalCombo->currentIndex() );

	 settings.setValue( "CP2K/NMultiGrid", ui.nMGridSpin->value() );
	 settings.setValue( "CP2K/CutOff", ui.cutOffSpin->value() );

	 // DFTB tab
	 settings.setValue( "CP2K/SCCType", ui.sccTypeCombo->currentIndex() );
	 settings.setValue( "CP2K/Dispersion", ui.dispersionCheck->checkState() );

	 // SE tab
	 settings.setValue( "CP2K/SEMethod", ui.seMethodCombo->currentIndex() );

	 // QMMM tab
	 settings.setValue( "CP2K/Ecouple", ui.qmmmElectroCoupleCombo->currentIndex() );

	 // MD tab
	 settings.setValue( "CP2k/Ensemble", ui.ensembleCombo->currentIndex() );
	 
	 // Common
	 settings.setValue("CP2K/SavePath", m_savePath);

  }


  void Cp2kInputDialog::readSettings(QSettings &settings)
  {
	 // Basic tab
	 QVariant defaultProjectName("myProject");
     m_projectName = settings.value("CP2K/ProjectName", defaultProjectName ).toString();
     ui.projectNameLine->setText(m_projectName);
	  
	 ui.runTypeCombo->setCurrentIndex(settings.value("CP2K/RunType", 0).toInt()); 
	 setRunType( settings.value( "CP2K/RunType", 0).toInt() );

	 ui.mmRadioButton->setChecked( settings.value("CP2K/MMRadio", true).toBool() );
	 ui.qmRadioButton->setChecked( settings.value("CP2K/QMRadio", false).toBool() );
	 ui.qmmmRadioButton->setChecked( settings.value("CP2K/QMMMRadio", false).toBool() );

	 if( settings.value("CP2K/MMRadio", true).toBool() ) mmRadioChecked();
	 else if( settings.value("CP2K/QMRadio", true).toBool() ) qmRadioChecked();
	 else qmmmRadioChecked();

	 // MM tab
	 ui.emaxSplineDoubleSpin->setValue(settings.value("CP2K/EmaxSpline",0.5 ).toDouble() );
	 setEmaxSpline( settings.value("CP2K/EmaxSpline",0.5 ).toDouble() );

	 ui.ewaldTypeCombo->setCurrentIndex( settings.value("CP2K/EwaldType", 0 ).toInt() );
	 setEwaldType( settings.value("CP2K/EwaldType", 0 ).toInt() );

	 // QM tab
	 ui.qmMethodCombo->setCurrentIndex(settings.value("CP2K/QmMethod", 0).toInt() );
	 setQmMethod( settings.value("CP2K/QmMethod", 0).toInt() );

	 ui.scfGuessCombo->setCurrentIndex(settings.value("CP2K/SCFGuess", 0).toInt() );
	 setSCFGuess(settings.value("CP2K/SCFGuess", 0).toInt());

	 ui.chargeSpin->setValue(settings.value("CP2K/Charge", 0 ).toInt() );
	 setCharge( settings.value("CP2K/Charge", 0 ).toInt() );

	 ui.multiplicitySpin->setValue(settings.value("CP2K/Multiplicity", 1).toInt() );
	 setMultiplicity(settings.value("CP2K/Multiplicity", 1).toInt());

	 ui.maxSCFSpin->setValue(settings.value("CP2K/MaxSCF", 50).toInt() );
	 setMaxSCF( settings.value("CP2K/MaxSCF", 50).toInt() );

	 // DFT tab
	 ui.basisSetCombo->setCurrentIndex( settings.value("CP2K/BasisSet", 0 ).toInt() );
	 setBasisSet( settings.value("CP2K/BasisSet", 0 ).toInt() );

	 ui.functionalCombo->setCurrentIndex( settings.value("CP2K/Functional", 0 ).toInt() );
	 setFunctional( settings.value("CP2K/Functional", 0 ).toInt() );

	 ui.nMGridSpin->setValue( settings.value("CP2K/NMultiGrid", 4 ).toInt() );
	 setNMultiGrid( settings.value("CP2K/NMultiGrid", 4 ).toInt() );

	 ui.cutOffSpin->setValue( settings.value("CP2K/CutOff", 30 ).toInt() );
	 setCutOff(settings.value("CP2K/CutOff", 30 ).toInt());

	 // DFTB tab
	 ui.sccTypeCombo->setCurrentIndex( settings.value("CP2K/SCCType").toInt() );
	 setDftbSCCType( settings.value("CP2K/SCCType").toInt() );

	 ui.dispersionCheck->setCheckState( (Qt::CheckState)settings.value("CP2K/Dispersion", 0).toInt() );
	 setDftbDisp( settings.value("CP2K/Dispersion", 0).toInt() );

	 // SE tab
	 ui.seMethodCombo->setCurrentIndex( settings.value("CP2K/SEMethod", 0).toInt() );
	 setSEMethod( settings.value("CP2K/SEMethod", 0).toInt() );

	 // QMMM tab
	 ui.qmmmElectroCoupleCombo->setCurrentIndex( settings.value("CP2k/Ecouple", 0).toInt() );
	 setQMMMElectrostaticCoupling( settings.value("CP2k/Ecouple", 0).toInt() );

	 // MD tab
	 ui.ensembleCombo->setCurrentIndex( settings.value("CP2K/Ensemble", 0).toInt() );
	 setEnsemble( settings.value("CP2K/Ensemble", 0).toInt() );

	 // Common
	 m_savePath = settings.value("CP2K/SavePath").toString();

  }


  QString Cp2kInputDialog::generateInputDeck()
  {
	  QString buffer;
	  QTextStream mol(&buffer);

	  QString molBaseName;

	  if( m_molecule != NULL )
	  {
		  	QFileInfo molNameInfo( m_molecule->fileName() );
	        molBaseName =  molNameInfo.baseName();
	  }

	  mol << "&GLOBAL\n";

     	  mol << " PROJECT " << m_projectName << "\n";
	      mol << " RUN_TYPE " << m_runType << "\n";

	  mol << "&END GLOBAL\n\n";


	  mol << "&FORCE_EVAL\n";

	  // MM or QMMM
	    if( m_mmRadioChecked || m_qmmmRadioChecked ) 
		{
			if(m_mmRadioChecked)
				mol << " METHOD FIST\n";
			else if(m_qmmmRadioChecked)
				mol << " METHOD QMMM\n";

		    mol << " &MM\n";

				mol << "  &FORCEFIELD\n";

				if( molBaseName != NULL )
				{ mol << "   parm_file_name  " << molBaseName << ".pot\n";}
				else { mol << "   parm_file_name untitled.pot\n";}
				mol << "   parmtype CHM\n";

				mol << "   IGNORE_MISSING_CRITICAL_PARAMS T\n";

				mol << "   &SPLINE\n";
				mol << "     EMAX_SPLINE " << m_emaxSpline << "\n";
				mol << "   &END SPLINE\n";

				mol << "  &END FORCEFIELD\n";

				mol << "  &POISSON\n";
				mol << "   &EWALD\n";
				mol << "     EWALD_TYPE " << m_ewaldType << "\n";
				mol << "   &END EWALD\n";
				mol << "  &END POISSON\n";

			mol << " &END MM\n";
		}
		// MM or QMMM ends.

		// QM or QMMM
		if( m_qmRadioChecked || m_qmmmRadioChecked ) 
		{
			mol << endl;

			if(m_qmRadioChecked)
				mol << " METHOD QS\n";

			// else if(m_qmmmRadioChecked) then use METHOD QMMM as the previous shown.

			// DFT
			if( m_qmMethod == "DFT" )
			{
			   mol << tr("# Please copy the GTH_BASIS_SETS and POTENTIAL files\n");
			   mol << tr("# which may have been in cp2k-x.x.x/data or cp2k-x.x.x/tests/QS directory\n");
			   mol << tr("# to the CP2K executable directory if you use these default files.\n");

			   mol << " &DFT\n";
			   mol << "   CHARGE " << m_charge << "\n";
			   mol << "   MULTIPLICITY " << m_multiplicity << "\n";

			   mol << "  BASIS_SET_FILE_NAME  GTH_BASIS_SETS\n";
			   mol << "  POTENTIAL_FILE_NAME  POTENTIAL\n";
			   mol << "  &MGRID\n";
			   mol << "   NGRIDS " << m_nMultiGrid << "\n";
			   mol << "   CUTOFF " << m_cutOff << "\n";
			   mol << "  &END MGRID\n";

			   mol << "  &XC\n";
               mol << "   &XC_FUNCTIONAL " << m_functional << "\n"; 
               mol << "   &END XC_FUNCTIONAL\n";
               mol << "  &END XC\n";

			   mol << " &END DFT\n";
			}

			else if (m_qmMethod == "DFTB")
			{
			   if( m_isScc )
			   {
			      mol << tr("# Before using DFTB-SCC, please copy the scc folder\n");
                  mol << tr("# which may have been in cp2k-x.x.x/data/DFTB or cp2k-x.x.x/tests/DFTB\n" );
			      mol << tr("# and its contents to ../(cp2k executable directory).\n");
				  mol << "#\n";
				  mol << tr("# If you use dispersion correction, please copy uff_table\n");
                  mol << tr("# which may have been in cp2k-x.x.x/data/DFTB or cp2k-x.x.x/tests/DFTB\n" );
			      mol << tr("# to ../(cp2k executable directory).\n");
			   }

			   else
			   {
				  mol << tr("# Before using DFTB-NONSCC, please copy the nonscc folder\n");
			      mol << tr("# which may have been in cp2k-x.x.x/data/DFTB or cp2k-x.x.x/tests/DFTB\n" );
			      mol << tr("# and its contents to ../(cp2k executable directory).\n");
				  mol << "#\n";
				  mol << tr("# If you use dispersion correction, please copy uff_table\n");
                  mol << tr("# which may have been in cp2k-x.x.x/data/DFTB or cp2k-x.x.x/tests/DFTB\n" );
			      mol << tr("# to ../(cp2k executable directory).\n");
			   }

			   mol << " &DFT\n";
			   mol << "   CHARGE " << m_charge << "\n";
			   mol << "   MULTIPLICITY " << m_multiplicity << "\n";
			   mol << "  &QS\n";
			   mol << "   METHOD DFTB\n";

			   mol << "   &DFTB\n";

			   if( m_isScc )
			   {  mol << "     SELF_CONSISTENT    T\n"; }
			   else
			   { mol << "     SELF_CONSISTENT    F\n";}

			   if( m_dispersion )
			   { mol << "     DISPERSION         T\n";}
			   else
			   { mol << "     DISPERSION         F\n";}

               mol << "     ORTHOGONAL_BASIS   F\n";  
               mol << "     DO_EWALD           F\n";
               mol << "     &PARAMETER\n";

			   if( m_isScc )
			   {
                   mol << "      PARAM_FILE_PATH  ../scc\n";
                   mol << "      PARAM_FILE_NAME  scc_parameter\n";
			   }
			   else
			   {
				   mol << "      PARAM_FILE_PATH  ../nonscc\n";
                   mol << "      PARAM_FILE_NAME  nonscc_parameter\n";
			   }

			   if( m_dispersion )
			   { mol << "      UFF_FORCE_FIELD  uff_table\n"; }

               mol << "     &END PARAMETER\n";
			   mol << "   &END DFTB\n";

			   mol << "  &END QS\n";

			   mol << "  &SCF\n";
			   mol << "   SCF_GUESS " << m_scfGuess << "\n";
			   mol << "   MAX_SCF " << m_maxSCF << "\n";
			   mol << "  &END SCF\n";

			   mol << " &END DFT\n";
			}

			else if (m_qmMethod == "SE")
			{
			   mol << " &DFT\n";
			   mol << "   CHARGE " << m_charge << "\n";
			   mol << "   MULTIPLICITY " << m_multiplicity << "\n";
			   mol << "  &QS\n";
			   mol << "   METHOD " << m_seMethod << "\n";
			   mol << "   &SE\n";
			   mol << "     ANALYTICAL_GRADIENTS F\n";
			   mol << "   &END SE\n";
			   mol << "  &END QS\n";

               mol << "  &SCF\n";
               mol << "    SCF_GUESS " << m_scfGuess << "\n";
			   mol << "    MAX_SCF  " << m_maxSCF << "\n";
               mol << "  &END SCF\n";

			   mol << " &END DFT\n";
			}

		}
		// QM or QMMM ends.

		// QMMM
		if( m_qmmmRadioChecked ) 
		{
		    mol << endl;
		    mol << " &QMMM\n";
			mol << "  ECOUPL " << m_qmmmElectrostaticCoupling << "\n";

			setSelAtoms();

			setLinkAtoms();

			if (!linkAtoms.empty())
			{
				for (int i = 0; i < linkAtoms.size(); i++)
				{
					mol << "  &LINK\n";
					mol << "   ALPHA 1.50\n";
					mol << "   LINK_TYPE IMOMM\n";
		            mol << "   MM_INDEX " << linkAtoms[i].mmidx << "\n";
					mol << "   QM_INDEX " << linkAtoms[i].qmidx << "\n";
					mol << "   QMMM_SCALE_FACTOR 0.0\n";
					mol << "   RADIUS 0.80\n";
					mol << "  &END LINK\n";
				}

			}

			if( !selAtoms.empty() && !selAtomsKind.empty() )
			{
				for (int i = 0; i < selAtomsKind.size(); i++)
				{
					mol << "  &QM_KIND " << selAtomsKind[i] << endl;
					mol << "   MM_INDEX";
					for (int j = 0; j < selAtoms.size(); j++)
					{
						if (selAtomsKind[i] == selAtoms[j].element)
						{
							mol << " " << selAtoms[j].idx;
						}
					}

					mol << endl;

					mol << "  &END QM_KIND\n";
				}
			}


		   mol << "  &CELL\n";
		   mol << unitCell();
		   mol << "  &END CELL\n";

		   mol << " &END QMMM\n\n";

       }

       // QMMM ends

		mol << " &SUBSYS\n";

		 // Define coord in input only for QM (Specify coord file for MM)
		if( (m_molecule != NULL) && m_qmRadioChecked )
		{
			  mol << "  &COORD\n";

			  QList<Atom *> atoms = m_molecule->atoms();

			  foreach (Atom *atom, atoms) 
			  {
			  mol << "    "
				  << qSetFieldWidth(3) << left << QString(OpenBabel::etab.GetSymbol(atom->atomicNumber()))
				  << qSetFieldWidth(15) << qSetRealNumberPrecision(5) << forcepoint << fixed << right 
			      << atom->pos()->x() << atom->pos()->y() << atom->pos()->z()
				  << qSetFieldWidth(0) << "\n";
			  }

			  mol << "  &END COORD\n\n";

		 }

		// Specify coord file for MM
		else if( m_mmRadioChecked || m_qmmmRadioChecked )
		{
			mol << "   &TOPOLOGY\n";

			if( molBaseName != NULL)
			{ 
				mol << "    CONN_FILE_NAME  " << molBaseName << ".psf\n";
				mol <<  "    CONNECTIVITY UPSF\n";
				mol << "    COORD_FILE_NAME  " << molBaseName << ".xyz\n";
				mol <<  "    COORDINATE XYZ\n";
			}

			else
			{
				mol << "    CONN_FILE_NAME untitled.psf\n";
				mol <<  "    CONNECTIVITY UPSF\n";
				mol << "    COORD_FILE_NAME untitled.xyz\n";
				mol <<  "    COORDINATE XYZ\n";
			}
	
			mol << "   &END TOPOLOGY\n";
		}


		// Atom kind for QM-DFT
		 if( m_molecule != NULL && (&atomKindMol != NULL) && (m_qmRadioChecked) && (m_qmMethod == "DFT") )
		 {
			 setAtomKindMol();

		     for( int i = 0 ; i < atomKindMol.size() ; i++ )
		     {
		         mol << "   &KIND " << atomKindMol[i] << "\n";

				 mol  << "    BASIS_SET " << m_basisSet << "\n";
				 mol  << "    POTENTIAL " 
					  << potentialName( atomKindMol[i] )
                      << "\n";

			     mol << "   &END KIND\n\n";
		      }
		 }
		 // Atom kind for QM-DFT ends.


		 // Atom kind for QMMM-DFT
		 if(  m_molecule != NULL && m_qmmmRadioChecked && (m_qmMethod == "DFT") && !selAtomsKind.empty() )
		 {
			 for (int i = 0; i < selAtomsKind.size(); i++)
			 {
				 mol << "   &KIND " << selAtomsKind[i] << "\n";

				 mol << "    BASIS_SET " << m_basisSet << "\n";
				 mol << "    POTENTIAL "
					 << potentialName(selAtomsKind[i])
					 << "\n";

				 mol << "   &END KIND\n";
			 }
		 }
		 // Atom kind for QMMM-DFT ends.

		mol << "  &CELL\n";
	    mol << unitCell();
	   	mol << "  &END CELL\n";

        /*
		if( Constraints != NULL )
		{
			mol << "   &COLVAR\n";

			for( int i = 0 ; i < Constraints->rowCount() ; i++ )
			{
				if( Constraints->constraints().GetConstraintType( i ) == OBFF_CONST_IGNORE )
				{
					mol << "#   atom " << Constraints->constraints().GetConstraintAtomA( i ) << " is ignored\n";
				}
			}

			mol << "   &END COLVAR\n";
		}
		*/


		mol << " &END SUBSYS\n";

	  mol << "&END FORCE_EVAL\n";

	  if( m_runType == "MD" )
	  {
		  mol << "&MOTION\n";
		  mol << " &MD\n";

          mol << "   ENSEMBLE " << m_ensemble << "\n";
          mol << "   STEPS 10\n";
          mol << "   TIMESTEP 0.5\n";
          mol << "   TEMPERATURE 298\n";

		  mol << " &END MD\n";
		  mol << "&END MOTION\n";
	  }

	  return buffer;

  }

  bool Cp2kInputDialog::qstringEqual(const QString& left, const QString& right)
  {
	  if( left == right ) return true;
	  else return false;
  }


  void Cp2kInputDialog::setSelAtoms()
  {
	  if( m_molecule == NULL ) return;

	  selAtoms.clear();

	  selatom tempsa;

	 // Find selected atoms
	 GLWidget *widget = GLWidget::current();// a pointer to the current GLWidget

	 if( widget->selectedPrimitives().size() )
     {
		 foreach(Primitive *p, widget->selectedPrimitives())
         {
              if (p->type() == Primitive::AtomType)
               {
                  Atom *a = static_cast<Atom *>(p);

				  tempsa.element = QString(etab.GetSymbol(a->atomicNumber()));
				  tempsa.uid = a->OBAtom().GetId() + 1;
				  tempsa.idx = m_molecule->OBMol().GetAtomById( a->OBAtom().GetId() )->GetIdx();

                  selAtoms.push_back(tempsa);
				}
  
          } 
	  }

	 selAtomsKind.clear();

	  for( int i=0; i< selAtoms.size(); i++)
		  selAtomsKind.push_back(selAtoms[i].element);

	  std::sort( selAtomsKind.begin(), selAtomsKind.end() );
	  std::vector<QString>::iterator new_end = std::unique( selAtomsKind.begin(), selAtomsKind.end(), qstringEqual );
	  selAtomsKind.erase( new_end, selAtomsKind.end());

  }

  QString Cp2kInputDialog::unitCell()
  {
	  QString buff;
	  QTextStream ts(&buff);

	  if (m_molecule != NULL)
	  {
		  OpenBabel::OBUnitCell* cell;
		  cell = m_molecule->OBUnitCell();

		  if (cell != NULL)
		  {
			  // lengths of cell vectors
			  ts << "    " << qSetFieldWidth(31) << left << "ABC" << qSetFieldWidth(0)
				  << qSetFieldWidth(15) << left << qSetRealNumberPrecision(5) << forcepoint << fixed
				  << cell->GetA() << cell->GetB() << cell->GetC()
				  << qSetFieldWidth(0) << "\n";

			  // angles between cell vectors
			  ts << "    " << qSetFieldWidth(20) << left << "ALPHA_BETA_GAMMA"
				  << qSetFieldWidth(15) << left << qSetRealNumberPrecision(5) << forcepoint << fixed
				  << cell->GetAlpha() << cell->GetBeta() << cell->GetGamma()
				  << qSetFieldWidth(0) << "\n";
		  }

		  else
		  {
			  ts << tr("   ABC 50 50 50 # Provisional Cell Param. Please \"Add Unit Cell\" !\n");
		  }

	  }

	  return buff;
  }

  void Cp2kInputDialog::setLinkAtoms()
  {
	  if (m_molecule == NULL) return;
	  if (selAtoms.empty() ) return;

	  // collect mm uids
	  std::vector<unsigned long> mmuids;

	  QList<Atom *> atoms = m_molecule->atoms();

	  foreach(Atom *atom, atoms)
	  {
		  bool is_mm = true;

		  for (int i = 0; i < selAtoms.size(); i++) // scan qm atoms
		  {
			  if (selAtoms[i].uid == atom->OBAtom().GetId() + 1 )
				  is_mm = false;
		  }

		  if (is_mm) mmuids.push_back(atom->OBAtom().GetId() + 1);

	  }

	  if (mmuids.empty()) return;

	  linkAtoms.clear();

	  // find mm atoms connected to a qm atom and pack pairs into linkAtoms.
	  for (int i = 0; i < selAtoms.size(); i++) // scan qm atoms
	  {
		  unsigned long qm = selAtoms[i].uid - 1;

		  Atom* a = m_molecule->atomById(qm);

		  foreach(unsigned long u, a->neighbors())
		  {
			  for (int j = 0; j < mmuids.size(); j++)
			  {
				  if (u + 1 == mmuids[j])
				  {
					  linkatoms templa;
					  templa.mmuid = mmuids[j];
					  templa.qmuid = qm + 1;
					  templa.mmidx = m_molecule->OBMol().GetAtomById(u)->GetIdx();
					  templa.qmidx = m_molecule->OBMol().GetAtomById(qm)->GetIdx();

					  linkAtoms.push_back(templa);
				  }
			  }
		  }


	  }

  }


  void Cp2kInputDialog::updatePreviewText()
  {
      ui.previewText->setText(generateInputDeck());
  }

  void Cp2kInputDialog::resetClicked()
  {
	  // Confirmation need!

	  // Basic Tab
	  m_projectName = "myProject";
      ui.projectNameLine->setText( m_projectName );

	  // m_runType = "ENERGY";
	  ui.runTypeCombo->setCurrentIndex(0);
	  setRunType(0);

	  ui.mmRadioButton->setChecked(true);
	  ui.qmRadioButton->setChecked(false);
	  ui.qmmmRadioButton->setChecked(false);
	  mmRadioChecked();

	  // MM tab
	  m_emaxSpline = 0.5;
	  ui.emaxSplineDoubleSpin->setValue( m_emaxSpline );

	  // m_ewaldType = "NONE";
	  ui.ewaldTypeCombo->setCurrentIndex(0);
	  setEwaldType( 0 );

	  // QM tab
	  // m_qmMethod = "DFT";
	  ui.qmMethodCombo->setCurrentIndex(0);
	  setQmMethod(0);

	  ui.scfGuessCombo->setCurrentIndex(0);
	  setSCFGuess(0);

	  m_charge = 0;
	  ui.chargeSpin->setValue( m_charge );

	  m_multiplicity = 1;
	  ui.multiplicitySpin->setValue( m_multiplicity );

	  m_maxSCF = 50;
	  setMaxSCF( m_maxSCF );

	  // DFT tab
	  // m_basisSet = "SZV-GTH";
	  ui.basisSetCombo->setCurrentIndex(0);
	  setBasisSet(0);

	  // m_functional = "BLYP";
	  ui.functionalCombo->setCurrentIndex(0);
	  setFunctional(0);

	  m_nMultiGrid = 4;
	  ui.nMGridSpin->setValue( m_nMultiGrid );

	  m_cutOff = 30;
	  ui.cutOffSpin->setValue( m_cutOff) ;

	  // DFTB tab
	  ui.sccTypeCombo->setCurrentIndex(0);
	  setDftbSCCType(0);

	  ui.dispersionCheck->setCheckState((Qt::CheckState)0);
	  setDftbDisp(0);

	  // SE tab
	  ui.seMethodCombo->setCurrentIndex(0);
	  setSEMethod(0);

	  // QMMM tab
	  ui.qmmmElectroCoupleCombo->setCurrentIndex(0);
	  setQMMMElectrostaticCoupling(0);

	  // MD tab
	  ui.ensembleCombo->setCurrentIndex(0);
	  setEnsemble(0);

  }

  void Cp2kInputDialog::generateClicked()
  {

	 saveInputFile( ui.previewText->toPlainText(),
                          tr("Cp2k Input File"), QString("inp"));
 
  }

  void Cp2kInputDialog::closeClicked()
  {
	  QSettings settings;
      writeSettings(settings);

	  close();
  }

  void Cp2kInputDialog::setProjectName()
  {
	m_projectName = ui.projectNameLine->displayText();
    updatePreviewText();
  
  }

  void Cp2kInputDialog::setRunType(int n)
  {
	  switch(n)
	  {
	    case 0:
		default:
			m_runType = "ENERGY";
			ui.mdTab->setEnabled(false);
			break;
	    case 1:
			m_runType = "ENERGY_FORCE";
			ui.mdTab->setEnabled(false);
			break;
		case 2:
			m_runType = "MD";
			ui.mdTab->setEnabled(true);
			break;
		case 3:
			m_runType = "GEO_OPT";
			ui.mdTab->setEnabled(false);
			break;
	  }

	  updatePreviewText();

  }

  void Cp2kInputDialog::setEwaldType(int n)
  {
	  switch (n)
	  {
	  case 0:
	  default:
		  m_ewaldType = "NONE";
		  break;
	  case 1:
		  m_ewaldType = "EWALD";
		  break;
	  case 2:
		  m_ewaldType = "PME";
		  break;
	  case 3:
		  m_ewaldType = "SPME";
		  break;
	  }

	  updatePreviewText();
  }

  void Cp2kInputDialog::setQmMethod(int n)
  {
	  switch(n)
	  {
	    case 0:
		default:
			m_qmMethod = "DFT";

			if( m_qmRadioChecked || m_qmmmRadioChecked )
			{
	          ui.qmTab->setEnabled(true);
	          ui.dftTab->setEnabled(true);
	          ui.dftbTab->setEnabled(false);
	          ui.seTab->setEnabled(false);
			}

			break;
	    case 1:
			m_qmMethod = "DFTB";

			if( m_qmRadioChecked || m_qmmmRadioChecked )
			{
	           ui.qmTab->setEnabled(true);
	           ui.dftTab->setEnabled(false);
	           ui.dftbTab->setEnabled(true);
	           ui.seTab->setEnabled(false);
			}

			break;
		case 2:
			m_qmMethod = "SE";

			if( m_qmRadioChecked || m_qmmmRadioChecked )
			{
	          ui.qmTab->setEnabled(true);
	          ui.dftTab->setEnabled(false);
	          ui.dftbTab->setEnabled(false);
	          ui.seTab->setEnabled(true);
			}

			break;

	  }

	  updatePreviewText();
  }

  void Cp2kInputDialog::setSCFGuess( int n )
  {
	  switch(n)
	  {
	    case 0:
		default:
			m_scfGuess = "ATOMIC";
			break;
		case 1:
			m_scfGuess = "RESTART";
			break;
		case 2:
			m_scfGuess = "RANDOM";
			break;
		case 3:
			m_scfGuess = "CORE";
			break;
		case 4:
			m_scfGuess = "DENSITIES";
			break;
		case 5:
			m_scfGuess = "HISTORY_RESTART";
			break;
		case 6:
			m_scfGuess = "MOPAC";
			break;
		case 7:
			m_scfGuess = "SPARSE";
			break;
		case 8:
			m_scfGuess = "NONE";
			break;
	  }

	  updatePreviewText();
  }

  void Cp2kInputDialog::setBasisSet(int n)
  {
	  switch(n)
	  {
	    case 0:
		default:
		    m_basisSet = "SZV-GTH";
		    break;
		case 1:
			m_basisSet = "DZV-GTH";
		    break;
		case 2:
			m_basisSet = "DZVP-GTH";
            break;
		case 3:
			m_basisSet = "TZVP-GTH";
			break;
		case 4:
			m_basisSet = "TZV2P-GTH";
			break;
	  }

	  updatePreviewText();
  }

  void Cp2kInputDialog::setFunctional(int n)
  {
	  switch(n)
	  {
	     case 0:
		 default:
			 m_functional = "BLYP";
			 break;
		 case 1:
		     m_functional = "BP";
			 break;
		 case 2:
			 m_functional = "HCTH120";
			 break;
		 case 3:
			 m_functional = "PADE";
			 break;
		 case 4:
			 m_functional = "PBE";
			 break;

	  }

     updatePreviewText();

  }

  void Cp2kInputDialog::setDftbSCCType(int n)
  {
	  switch(n)
	  {
	    case 0:
		default:
		  m_isScc = true;
		  break;
		case 1:
		  m_isScc = false;
		  break;
	  }

	  updatePreviewText();
  }


  void Cp2kInputDialog::setDftbDisp(int n)
  {
	  switch(n)
	  {
	    case 0:
		default:
	      m_dispersion = false;
		  break;
		case 1:
		case 2:
		  m_dispersion = true;
		  break;
	  }

	 updatePreviewText();
  }

  void Cp2kInputDialog::setSEMethod(int n)
  {
	  switch(n)
	  {
	     case 0:
		 default:
			 m_seMethod = "PM6";
			 break;
		 case 1:
		     m_seMethod = "MNDO";
			 break;
		 case 2:
			 m_seMethod = "MNDOD";
			 break;
		 case 3:
			 m_seMethod = "AM1";
			 break;
		 case 4:
			 m_seMethod = "PM3";
			 break;
		 case 5:
			 m_seMethod = "PDG";
			 break;
		 case 6:
			 m_seMethod = "RM1";
			 break;
		 case 7:
			 m_seMethod = "PNNL";
			 break;
	  }

     updatePreviewText();

  }

  void Cp2kInputDialog::setQMMMElectrostaticCoupling(int n)
  {
	  switch(n)
	  {
	    case 0:
		default:
			m_qmmmElectrostaticCoupling = "NONE";
			break;
		case 1:
            m_qmmmElectrostaticCoupling = "COULOMB";
			break;
		case 2:
            m_qmmmElectrostaticCoupling = "GAUSS";
			break;
		case 3:
			m_qmmmElectrostaticCoupling = "S-WAVE";
            break;
	  }

	  updatePreviewText();
  }

   void Cp2kInputDialog::setEnsemble(int n)
   {
	   switch(n)
	   {
	     case 0:
	     default:
		     m_ensemble = "NVE";
		     break;
		 case 1:
			 m_ensemble = "NVT";
			 break;
	   }

	    updatePreviewText();
   }
  

  // Extracts atom kinds from m_molecule.
  void Cp2kInputDialog::setAtomKindMol()
  {
	  atomKindMol.clear();

	  if( m_molecule != NULL )
	  {
		  QList<Atom *> atoms = m_molecule->atoms();

		  std::vector<int> v0;

		  foreach (Atom *atom, atoms) 
		  { v0.push_back( atom->atomicNumber() ); }

		  // delete duplicate atoms.
		  sort(v0.begin(), v0.end());

		  std::vector<int>::iterator it;
		  it = unique(v0.begin(), v0.end());
		  v0.erase(it, v0.end());

		  for( int i=0; i<v0.size(); i++)
		  { atomKindMol.push_back(OpenBabel::etab.GetSymbol(v0[i]));}
	  }

	  return;
  }


  QString Cp2kInputDialog::potentialName( QString atomType )
  {

	  QString potName("GTH-");
	  potName.append( m_functional );
	  potName.append( "-q" );

	  int iAtomicNum = etab.GetAtomicNum( atomType.toStdString().c_str() );
      int iValenceElec = 0;

     if( iAtomicNum >= 1 && iAtomicNum <= 4 ) // From H to Be
	    { iValenceElec = iAtomicNum; }

	 else if( iAtomicNum >= 5 && iAtomicNum <= 12 ) // From B to Mg
	    { iValenceElec = iAtomicNum - 2;}

	 else if( iAtomicNum >= 13 && iAtomicNum <= 28 ) // From Al to Ni
	    { iValenceElec = iAtomicNum - 10;}

	 else if( iAtomicNum >= 29 && iAtomicNum <= 31 ) // Cu, Zn, Ga
        { iValenceElec = iAtomicNum - 18;}

	 else if( iAtomicNum >= 32 && iAtomicNum <= 36 ) // Ge, As, Br, Kr
		 { iValenceElec = iAtomicNum - 28;}

	  if( iValenceElec == 0 )
	  {
		  potName = tr("# type pseudopotential functional name !");
		  return potName;
	  }

	  char cValenceElec[4];
	  sprintf( cValenceElec, "%d", iValenceElec );

      potName.append( cValenceElec );

	  return potName;
  }


  void Cp2kInputDialog::mmRadioChecked()
  {
	  m_mmRadioChecked = true;
	  m_qmRadioChecked = false;
	  m_qmmmRadioChecked = false;

	  /*
	  ui.qmMethodLabel->hide();
	  ui.qmMethodCombo->hide();
	  */

	  ui.mmTab->setEnabled(true);
	  ui.qmTab->setEnabled(false);
	  ui.dftTab->setEnabled(false);
	  ui.dftbTab->setEnabled(false);
	  ui.seTab->setEnabled(false);
	  ui.qmmmTab->setEnabled(false);

	  updatePreviewText();
  }

  void Cp2kInputDialog::qmRadioChecked()
  {
	  m_mmRadioChecked = false;
	  m_qmRadioChecked = true;
	  m_qmmmRadioChecked = false;

	  /*
	  ui.qmMethodLabel->show();
	  ui.qmMethodCombo->show();
	  */

	  ui.mmTab->setEnabled(false);

	  if(m_qmMethod == "DFT")
	  {
	    ui.qmTab->setEnabled(true);
	    ui.dftTab->setEnabled(true);
	    ui.dftbTab->setEnabled(false);
	    ui.seTab->setEnabled(false);
	  }

	  else if( m_qmMethod == "DFTB" )
	  {
		ui.qmTab->setEnabled(true);
	    ui.dftTab->setEnabled(false);
	    ui.dftbTab->setEnabled(true);
	    ui.seTab->setEnabled(false);
	  }

	  else if(m_qmMethod == "SE")
	  {
		ui.qmTab->setEnabled(true);
	    ui.dftTab->setEnabled(false);
	    ui.dftbTab->setEnabled(false);
	    ui.seTab->setEnabled(true);
	  }

	  ui.qmmmTab->setEnabled(false);

	  updatePreviewText();
  }

  void Cp2kInputDialog::qmmmRadioChecked()
  {
	  m_mmRadioChecked = false;
	  m_qmRadioChecked = false;
	  m_qmmmRadioChecked = true;

	  ui.mmTab->setEnabled(true);

	  /*
	  ui.qmMethodLabel->show();
	  ui.qmMethodCombo->show();
	  */

	  if(m_qmMethod == "DFT")
	  {
	    ui.qmTab->setEnabled(true);
	    ui.dftTab->setEnabled(true);
	    ui.dftbTab->setEnabled(false);
	    ui.seTab->setEnabled(false);
	  }

	  else if( m_qmMethod == "DFTB" )
	  {
		ui.qmTab->setEnabled(true);
	    ui.dftTab->setEnabled(false);
	    ui.dftbTab->setEnabled(true);
	    ui.seTab->setEnabled(false);
	  }

	  else if(m_qmMethod == "SE")
	  {
		ui.qmTab->setEnabled(true);
	    ui.dftTab->setEnabled(false);
	    ui.dftbTab->setEnabled(false);
	    ui.seTab->setEnabled(true);
	  }

	   ui.qmmmTab->setEnabled(true);

	  updatePreviewText();
  }



  QString Cp2kInputDialog::saveInputFile(QString inputDeck, QString fileType, QString ext)
  {
    // Try to set default save path for dialog using the next sequence:
    //  1) directory of current file (if any);
    //  2) directory where previous deck was saved;
    //  3) $HOME
    QFileInfo defaultFile;
    QString defaultPath;

    if( m_molecule != NULL )
	{
       defaultFile.setFile(m_molecule->fileName());
       defaultPath = defaultFile.canonicalPath();
	}

	else
	{
		defaultFile.setFile( tr("untitled.inp") );
		defaultPath = "";
	}

    if(m_savePath == "") 
	{
      if (defaultPath.isEmpty())
        defaultPath = QDir::homePath();
    } 
	
	else  defaultPath = m_savePath;

    QString defaultFileName = defaultPath + '/' + defaultFile.baseName() + "." + ext;

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save CP2K Input"),
                                                      defaultFileName, fileType + " (*." + ext + ")");
	m_savePath = defaultPath; // set current path to save path.

    if(fileName == "")
      return fileName;

    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text)) 
	   return QString();

    file.write(inputDeck.toLocal8Bit()); // prevent troubles in Windows
    file.close(); // flush buffer!
    m_savePath = QFileInfo(file).absolutePath();

    return fileName;
  }

  void Cp2kInputDialog::setEmaxSpline( double value )
  {
	  m_emaxSpline = value;
	  updatePreviewText();
  }

  void Cp2kInputDialog::setCharge( int n )
  {
	  m_charge = n;
	  updatePreviewText();
  }

  void Cp2kInputDialog::setMultiplicity( int n )
  {
	  m_multiplicity = n;
	  updatePreviewText();
  }

  void Cp2kInputDialog::setMaxSCF( int n )
  {
	  m_maxSCF = n;
	  updatePreviewText();
  }

  void Cp2kInputDialog::setNMultiGrid( int n )
  {
	  m_nMultiGrid = n;
	  updatePreviewText();
  }

   void Cp2kInputDialog::setCutOff( int n )
  {
	  m_cutOff = n;
	  updatePreviewText();
  }

}

