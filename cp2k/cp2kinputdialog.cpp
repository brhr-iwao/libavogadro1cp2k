/**********************************************************************
  Cp2kInpuDialog - CP2k Input Dialog

  Copyright (C) Aoyama Iwao

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.cc/>

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

	  connect(ui.projectNameLine, SIGNAL(editingFinished()),this, SLOT(setProjectName()) );
	  connect(ui.runTypeCombo, SIGNAL(currentIndexChanged(int)),this, SLOT(setRunType(int)));
	  connect(ui.closeButton, SIGNAL(clicked()), this, SLOT(close()) );
	  connect(ui.generateButton, SIGNAL(clicked()), this, SLOT(generateClicked()) );

	  connect(ui.mmRadioButton, SIGNAL(clicked()), this, SLOT(mmRadioChecked()) );
	  connect(ui.qmRadioButton, SIGNAL(clicked()), this, SLOT(qmRadioChecked()) );

	  m_runType = "ENERGY";

	  ui.mmRadioButton->setChecked(false);
	  ui.qmRadioButton->setChecked(true);
	  m_mmRadioChecked = false;
	  m_qmRadioChecked = true;

	  updatePreviewText();

  }

  Cp2kInputDialog::~Cp2kInputDialog()
  {

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


  void Cp2kInputDialog::readSettings(QSettings&)
  {

  }

  void Cp2kInputDialog::writeSettings(QSettings&) const
  {

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

	    if( m_mmRadioChecked ) // MM
		{
			mol << " METHOD FIST\n";
		    mol << " &MM\n";

				mol << "  &FORCEFIELD\n";

				if( molBaseName != NULL )
				{ mol << "   parm_file_name  " << molBaseName << ".pot\n";}
				else { mol << "   parm_file_name untitled.pot\n";}
				mol << "   parmtype CHM\n";

				mol << "   IGNORE_MISSING_CRITICAL_PARAMS T\n";

				mol << "  &END FORCEFIELD\n";

				mol << "  &POISSON\n";
				mol << "   &EWALD\n";
				mol << "     EWALD_TYPE NONE\n";
				mol << "   &END EWALD\n";
				mol << "  &END POISSON\n";

			mol << " &END MM\n";
		}

		else if( m_qmRadioChecked ) // QM
		{
			mol << " METHOD QUICKSTEP\n";

		}

		mol << " &SUBSYS\n";

		if( m_qmRadioChecked && (m_molecule != NULL)) // QM
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


		 setAtomKind();

		 if( m_molecule != NULL && (&atomKind != NULL) && (m_qmRadioChecked) )
		 {
			 int i;
		     for( i = 0 ; i < atomKind.size() ; i++ )
		     {
		         mol << "   &KIND " << atomKind[i] << "\n";

			     mol  << "    BASIS_SET\n"
					  << "    POTENTIAL\n";

			     mol << "   &END KIND\n\n";
		      }
		 }

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
				  mol << tr("# Please Add Unit Cell or type cell parameters !\n");
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

  }

  void Cp2kInputDialog::generateClicked()
  {

	 saveInputFile( ui.previewText->toPlainText(),
                          tr("Cp2k Input File"), QString("inp"));

  
  }

  void Cp2kInputDialog::setProjectName()
  {
	m_projectName = ui.projectNameLine->text();
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


  void Cp2kInputDialog::mmRadioChecked()
  {
	  m_mmRadioChecked = true;
	  m_qmRadioChecked = false;

	  ui.mmTab->show();
	  ui.qmTab->hide();

	  updatePreviewText();
  }

  void Cp2kInputDialog::qmRadioChecked()
  {
	  m_mmRadioChecked = false;
	  m_qmRadioChecked = true;

	  ui.mmTab->hide();
	  ui.qmTab->show();

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

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Input Deck"),
                                                      defaultFileName, fileType + " (*." + ext + ")");

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
}

