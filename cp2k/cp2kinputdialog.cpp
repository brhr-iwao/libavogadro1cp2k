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
                QDialog(parent, f), m_molecule(0),
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

	  // MM tab
	  connect(ui.emaxSplineDoubleSpin, SIGNAL(valueChanged(double)), this, SLOT(setEmaxSpline(double)));
	  connect(ui.ewaldTypeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(setEwaldType(int)));

	  // QM tab
	  connect(ui.qmMethodCombo,SIGNAL(currentIndexChanged(int)),this, SLOT(setQmMethod(int)) );
	  connect(ui.chargeSpin, SIGNAL(valueChanged(int)), this, SLOT(setCharge(int)));
	  connect(ui.multiplicitySpin, SIGNAL(valueChanged(int)), this, SLOT(setMultiplicity(int)));

	  // DFT tab
	  connect(ui.basisSetCombo,SIGNAL(currentIndexChanged(int)),this, SLOT(setBasisSet(int)) );
	  connect(ui.functionalCombo,SIGNAL(currentIndexChanged(int)),this, SLOT(setFunctional(int)) );
	  connect(ui.nMGridSpin, SIGNAL(valueChanged(int)), this, SLOT(setNMultiGrid(int)));
	  connect(ui.cutOffSpin, SIGNAL(valueChanged(int)), this, SLOT(setCutOff(int)));

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

	  // DFT tab
	  ui.nMGridSpin->setRange(1,10);
	  ui.nMGridSpin->setSingleStep(1);

	  ui.cutOffSpin->setRange(1, 1000);
	  ui.cutOffSpin->setSingleStep(1);

	  QSettings settings;
      readSettings(settings);

	  updatePreviewText();

  }

  Cp2kInputDialog::~Cp2kInputDialog()
  {
	  // the following codes don't work.
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

  void Cp2kInputDialog::setModel(ConstraintsModel *model)
  {
	  if (m_constraints) disconnect(m_constraints, 0, this, 0);
      m_constraints = model;
  }

  void Cp2kInputDialog::writeSettings(QSettings &settings) const
  {
    // Basic tab
	 settings.setValue("CP2K/ProjectName", ui.projectNameLine->displayText() );
     settings.setValue("CP2K/RunType", ui.runTypeCombo->currentIndex());
	 settings.setValue("CP2K/mmRadio", m_mmRadioChecked );
	 settings.setValue("CP2K/qmRadio", m_qmRadioChecked );

	 // MM tab
     settings.setValue("CP2K/EmaxSpline", ui.emaxSplineDoubleSpin->value() );
	 settings.setValue("CP2K/EwaldType", ui.ewaldTypeCombo->currentIndex() );

	 // QM tab
	 settings.setValue("CP2K/QmMethod", ui.qmMethodCombo->currentIndex() );
	 settings.setValue("CP2K/Charge", ui.chargeSpin->value() );
	 settings.setValue("CP2K/Multiplicity", ui.multiplicitySpin->value() );

	 // DFT tab
	 settings.setValue( "CP2K/BasisSet", ui.basisSetCombo->currentIndex() );
	 settings.setValue( "CP2K/Functional", ui.functionalCombo->currentIndex() );

	 settings.setValue( "CP2K/NMultiGrid", ui.nMGridSpin->value() );
	 settings.setValue( "CP2K/CutOff", ui.cutOffSpin->value() );

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

	 ui.mmRadioButton->setChecked( settings.value("CP2K/mmRadio", true).toBool() );
	 ui.qmRadioButton->setChecked( settings.value("CP2K/qmRadio", false).toBool() );

	 if( settings.value("CP2K/mmRadio", true).toBool() ) mmRadioChecked();
	 else qmRadioChecked();

	 // MM tab
	 ui.emaxSplineDoubleSpin->setValue(settings.value("CP2K/EmaxSpline",0.5 ).toDouble() );
	 setEmaxSpline( settings.value("CP2K/EmaxSpline",0.5 ).toDouble() );

	 ui.ewaldTypeCombo->setCurrentIndex( settings.value("CP2K/EwaldType", 0 ).toInt() );
	 setEwaldType( settings.value("CP2K/EwaldType", 0 ).toInt() );

	 // QM tab
	 ui.qmMethodCombo->setCurrentIndex(settings.value("CP2K/QmMethod", 0).toInt() );
	 setQmMethod( settings.value("CP2K/QmMethod", 0).toInt() );

	 ui.chargeSpin->setValue(settings.value("CP2K/Charge", 0 ).toInt() );
	 setCharge( settings.value("CP2K/Charge", 0 ).toInt() );

	 ui.multiplicitySpin->setValue(settings.value("CP2K/Multiplicity", 1).toInt() );
	 setMultiplicity(settings.value("CP2K/Multiplicity", 1).toInt());

	 // DFT tab
	 ui.basisSetCombo->setCurrentIndex( settings.value("CP2K/BasisSet", 0 ).toInt() );
	 setBasisSet( settings.value("CP2K/BasisSet", 0 ).toInt() );

	 ui.functionalCombo->setCurrentIndex( settings.value("CP2K/Functional", 0 ).toInt() );
	 setFunctional( settings.value("CP2K/Functional", 0 ).toInt() );

	 ui.nMGridSpin->setValue( settings.value("CP2K/NMultiGrid", 4 ).toInt() );
	 setNMultiGrid( settings.value("CP2K/NMultiGrid", 4 ).toInt() );

	 ui.cutOffSpin->setValue( settings.value("CP2K/CutOff", 30 ).toInt() );
	 setCutOff(settings.value("CP2K/CutOff", 30 ).toInt());

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

	  // MM
	    if( m_mmRadioChecked ) 
		{
			mol << " METHOD FIST\n";
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
		// MM ends.

		else if( m_qmRadioChecked ) // QM
		{

			// DFT
			if( m_qmMethod == "DFT" )
			{
			   mol << tr("# Please copy the GTH_BASIS_SETS and POTENTIAL files, \n");
			   mol << tr("# which may have been in cp2k-2.x.x/data or cp2k-2.x.x/tests/QS directory,\n");
			   mol << tr("# to the CP2K executable directory if you use these default files.\n");

			   mol << " METHOD QS\n";
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

			else if (m_qmMethod == "DFTB-SCC")
			{
			   mol << tr("# Please copy the scc folder,\n");
               mol << tr("# which may have been in cp2k-2.x.x/data/DFTB or cp2k-2.x.x/tests/DFTB,\n" );
			   mol << tr("# and its contents to ../(cp2k executable directory).\n");

			   mol << " &DFT\n";
			   mol << "   CHARGE " << m_charge << "\n";
			   mol << "   MULTIPLICITY " << m_multiplicity << "\n";
			   mol << "  &QS\n";
			   mol << "   METHOD DFTB\n";

			   mol << "   &DFTB\n";
               mol << "     SELF_CONSISTENT    T\n";
               mol << "     DISPERSION         F\n";
               mol << "     ORTHOGONAL_BASIS   F\n";  
               mol << "     DO_EWALD           F\n";
               mol << "     &PARAMETER\n";
               mol << "      PARAM_FILE_PATH  ../scc\n";
               mol << "      PARAM_FILE_NAME  scc_parameter\n";
               mol << "     &END PARAMETER\n";
			   mol << "   &END DFTB\n";

			   mol << "  &END QS\n";

			   mol << "  &SCF\n";
			   mol << "   SCF_GUESS CORE\n";
			   mol << "   MAX_SCF  20\n";
			   mol << "  &END SCF\n";

			   mol << " &END DFT\n";
			}

			else if (m_qmMethod == "DFTB-NONSCC")
			{
			   mol << tr("# Please copy the nonscc folder,\n");
			   mol << tr("# which may have been in cp2k-2.x.x/data/DFTB or cp2k-2.x.x/tests/DFTB,\n" );
			   mol << tr("# and its contents to ../(cp2k executable directory).\n");

			   mol << " &DFT\n";
			   mol << "   CHARGE " << m_charge << "\n";
			   mol << "   MULTIPLICITY " << m_multiplicity << "\n";
			   mol << "  &QS\n";
			   mol << "   METHOD DFTB\n";

			   mol << "   &DFTB\n";
               mol << "     SELF_CONSISTENT    T\n";
               mol << "     DISPERSION         F\n";
               mol << "     ORTHOGONAL_BASIS   F\n";  
               mol << "     DO_EWALD           F\n";
               mol << "     &PARAMETER\n";
               mol << "      PARAM_FILE_PATH  ../nonscc\n";
               mol << "      PARAM_FILE_NAME  nonscc_parameter\n";
               mol << "     &END PARAMETER\n";
			   mol << "   &END DFTB\n";

			   mol << "  &END QS\n";

			   mol << "  &SCF\n";
			   mol << "   SCF_GUESS CORE\n";
			   mol << "   MAX_SCF  20\n";
			   mol << "  &END SCF\n";

			   mol << " &END DFT\n";
			}

			else if (m_qmMethod == "SE-PM6")
			{
			   mol << " METHOD QS\n";
			   mol << " &DFT\n";
			   mol << "   CHARGE " << m_charge << "\n";
			   mol << "   MULTIPLICITY " << m_multiplicity << "\n";
			   mol << "  &QS\n";
			   mol << "   METHOD PM6\n";
			   mol << "   &SE\n";
			   mol << "     ANALYTICAL_GRADIENTS F\n";
			   mol << "   &END SE\n";
			   mol << "  &END QS\n";

               mol << "  &SCF\n";
               mol << "    SCF_GUESS ATOMIC\n";
               mol << "&END SCF\n";

			   mol << " &END DFT\n";
			}

		}

		mol << " &SUBSYS\n";

		 // Define coord in input only for QM (Specify coord file for MM)
		if( m_qmRadioChecked && (m_molecule != NULL))
		{
			  QList<Atom *> atoms = m_molecule->atoms();

			  mol << "  &COORD\n";

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
		else if( m_mmRadioChecked )
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


		// Atom kind for DFT
		 if( m_molecule != NULL && (&atomKind != NULL) && (m_qmRadioChecked) && (m_qmMethod == "DFT") )
		 {
			 setAtomKind();

			 int i;
		     for( i = 0 ; i < atomKind.size() ; i++ )
		     {
		         mol << "   &KIND " << atomKind[i] << "\n";

				 mol  << "    BASIS_SET " << m_basisSet << "\n";
				 mol  << "    POTENTIAL " 
					  << potentialName( atomKind[i] )
                      << "\n";

			     mol << "   &END KIND\n\n";
		      }
		 }
		 // Atom kind for DFT ends.

		 mol << "  &CELL\n";
  		 if(m_molecule != NULL)
		 {
			  OpenBabel::OBUnitCell* cell;
			  cell = m_molecule->OBUnitCell();

			  if(cell != NULL) 
			  {
				 // lengths of cell vectors
				  mol << "    " << qSetFieldWidth(31) << left << "ABC" << qSetFieldWidth(0) 
					 << qSetFieldWidth(15) << left <<  qSetRealNumberPrecision(5) << forcepoint << fixed
				     <<  cell->GetA() << cell->GetB() << cell->GetC() 
					 << qSetFieldWidth(0) << "\n";

				 // angles between cell vectors
				 mol << "    " << qSetFieldWidth(20) << left << "ALPHA_BETA_GAMMA"
					 << qSetFieldWidth(15) << left << qSetRealNumberPrecision(5) << forcepoint << fixed
				    <<  cell->GetAlpha() << cell->GetBeta() << cell->GetGamma() 
					 << qSetFieldWidth(0) << "\n";
			   }

			  else
			  {
				  mol << tr("   ABC 50 50 50 # Provisional Cell Param. Please \"Add Unit Cell\" !\n");
			  }

		 }
	   	mol << "  &END CELL\n";

		 mol << " &END SUBSYS\n";

	  mol << "&END FORCE_EVAL";


	  return buffer;

  }

  void Cp2kInputDialog::updatePreviewText()
  {
      ui.previewText->setText(generateInputDeck());
  }

  void Cp2kInputDialog::resetClicked()
  {
	  // Confirmation is need!

	  // Basic Tab
	  m_projectName = "myProject";
      ui.projectNameLine->setText( m_projectName );

	  // m_runType = "ENERGY";
	  ui.runTypeCombo->setCurrentIndex(0);
	  setRunType(0);

	  ui.mmRadioButton->setChecked(true);
	  ui.qmRadioButton->setChecked(false);
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

	  m_charge = 0;
	  ui.chargeSpin->setValue( m_charge );

	  m_multiplicity = 1;
	  ui.multiplicitySpin->setValue( m_multiplicity );

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
			m_runType = "ENERGY";
			break;
	    case 1:
			m_runType = "ENERGY_FORCE";
			break;
		case 2:
			m_runType = "MD";
			break;
		case 3:
			m_runType = "GEO_OPT";
			break;
		case 4:
			m_runType = "MC";
			break;
		default:
			m_runType = "ENERGY";
			break;

	  }

	  updatePreviewText();

  }

  void Cp2kInputDialog::setEwaldType(int n)
  {
	  switch (n)
	  {
	  case 0:
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
	  default:
		  m_ewaldType = "NONE";
		  break;
	  }

	  updatePreviewText();
  }

  void Cp2kInputDialog::setQmMethod(int n)
  {
	  switch(n)
	  {
	    case 0:
			m_qmMethod = "DFT";
			break;
	    case 1:
			m_qmMethod = "DFTB-SCC";
			break;
		case 2:
			m_qmMethod = "DFTB-NONSCC";
			break;
		case 3:
			m_qmMethod = "SE-PM6";
			break;
		default:
			m_qmMethod = "DFT";
			break;

	  }

	  updatePreviewText();
  }

  void Cp2kInputDialog::setBasisSet(int n)
  {
	  switch(n)
	  {
	    case 0:
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
		default:
			m_basisSet = "SZV-GTH";
			break;
	  }

	  updatePreviewText();
  }

  void Cp2kInputDialog::setFunctional(int n)
  {
	  switch(n)
	  {
	     case 0:
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
		 default:
			 m_functional = "BLYP";
			 break;
	  }

     updatePreviewText();

  }
  

  // Extracts atom kinds
  void Cp2kInputDialog::setAtomKind()
  {
	  atomKind.clear();

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

		  int i;
		  for( i=0; i<v0.size(); i++)
		  { atomKind.push_back(OpenBabel::etab.GetSymbol(v0[i]));}
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

	  /*
	  ui.qmMethodLabel->hide();
	  ui.qmMethodCombo->hide();
	  */

	  ui.mmTab->setEnabled(true);
	  ui.qmTab->setEnabled(false);
	  ui.dftTab->setEnabled(false);

	  updatePreviewText();
  }

  void Cp2kInputDialog::qmRadioChecked()
  {
	  m_mmRadioChecked = false;
	  m_qmRadioChecked = true;

	  /*
	  ui.qmMethodLabel->show();
	  ui.qmMethodCombo->show();
	  */

	  ui.mmTab->setEnabled(false);
	  ui.qmTab->setEnabled(true);
	  ui.dftTab->setEnabled(true);

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
	// m_savePath = defaultPath; // set current path to save path.

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

