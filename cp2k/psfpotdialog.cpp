/**********************************************************************
  PsfPotDialog - Charmm-style Psf/Pot Generator Dialog

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

#include "psfpotdialog.h"

#include <avogadro/molecule.h>
#include <avogadro/atom.h>
#include <avogadro/residue.h>

// openbabel headers were included in psfpotdialog.h
/*
#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/residue.h>
*/

#include <QString>
#include <QTextStream>
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>

using namespace OpenBabel;
using namespace std;

namespace Avogadro
{
  PsfPotDialog::PsfPotDialog(QWidget *parent, Qt::WindowFlags f ) : 
       QDialog(parent, f), m_molecule(0), m_savePath("")
  {
	  ui.setupUi(this);

	  connect(ui.ffCombo, SIGNAL(currentIndexChanged(int)),this, SLOT(setForceField(int)));

	  // buttons
	  connect(ui.closeButton, SIGNAL(clicked()), this, SLOT(closeClicked()) );
	  connect(ui.psfGenerateButton, SIGNAL(clicked()), this, SLOT( psfGenerateClicked() ) );
	  connect(ui.potGenerateButton, SIGNAL(clicked()), this, SLOT( potGenerateClicked() ) );
	  connect(ui.updatePreviewButton, SIGNAL(clicked()), this, SLOT( updatePreviewText() ) );

	  Updated = true;

	  QSettings settings;
	  readSettings( settings );

	  updatePreviewText();

  }


  PsfPotDialog::~PsfPotDialog()
  {
	  QSettings settings;
	  writeSettings( settings );
  }

  void PsfPotDialog::setMolecule(Molecule *molecule)
  {
	  // Disconnect the old molecule first...
	  if (m_molecule) disconnect(m_molecule, 0, this, 0);

	  Updated = false;
   
       m_molecule = molecule;
  }

  void PsfPotDialog::readSettings(QSettings &settings )
  {
	  ui.ffCombo->setCurrentIndex(settings.value("PsfPot/FF",0).toInt() );
	  setForceField(settings.value("PsfPot/FF",0).toInt());

	  m_savePath = settings.value("PsfPot/SavePath").toString();

  }

  void PsfPotDialog::writeSettings(QSettings &settings) const
  {
	  settings.setValue("PsfPot/FF", ui.ffCombo->currentIndex() );
	  settings.setValue("PsfPot/SavePath", m_savePath);
  }

 

  QString PsfPotDialog::psfPreviewPane()
  {
	  // GAFF case
	  if( m_forceField == "GAFF")
	       return psfGaff();

	  // UFF case
	  else if( m_forceField == "UFF")
	       return psfUff();

	  // Ghemical case
	  else if( m_forceField == "Ghemical")
	       return psfGhemical();

      // default
	  else return psfGaff();

  }

  QString PsfPotDialog::potPreviewPane()
  {
	  // GAFF case
	  if( m_forceField == "GAFF")
	       return potGaff();
	  
	  // UFF case
	  else if( m_forceField == "UFF")
		  return potUff();

	  // Ghemical case
	  else if( m_forceField == "Ghemical")
	       return potGhemical();

	  // default
	  else return potGaff();
  }

  void PsfPotDialog::closeClicked()
  {
	  QSettings settings;
      writeSettings(settings);

	  close();
  }

  void PsfPotDialog::updatePreviewText()
  {
      ui.psfPreviewText->setText(psfPreviewPane());
	  ui.potPreviewText->setText(potPreviewPane());

	  Updated = true;
  }

  void PsfPotDialog::psfGenerateClicked()
  {
	  if(!Updated)
	  {
		  QMessageBox msgBox;
          msgBox.setText(tr("The moleule has been modified."));
          msgBox.setInformativeText(tr("Will you update before generating file?"));
          msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
          msgBox.setDefaultButton(QMessageBox::Yes);
          int ret = msgBox.exec();

		  switch (ret) 
		  {
              case QMessageBox::Yes:
                 return;
              case QMessageBox::No:
                 break;
              default:
                 return;
          }
	  }

	 saveInputFile( ui.psfPreviewText->toPlainText(),
                          tr("Protein Structure File"), QString("psf"));
  
  }

  void PsfPotDialog::potGenerateClicked()
  {
	  if(!Updated)
	  {
		  QMessageBox msgBox;
          msgBox.setText(tr("The moleule has been modified."));
          msgBox.setInformativeText(tr("Will you update before generating file?"));
          msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
          msgBox.setDefaultButton(QMessageBox::Yes);
          int ret = msgBox.exec();

		  switch (ret) 
		  {
              case QMessageBox::Yes:
                 return;
              case QMessageBox::No:
                 break;
              default:
                 return;
          }
	  }

	 saveInputFile( ui.potPreviewText->toPlainText(),
                          tr("FF Parameter File"), QString("pot"));

  
  }

  QString PsfPotDialog::saveInputFile(QString inputDeck, QString fileType, QString ext)
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
		defaultFile.setFile("untitled");
		defaultPath = "";
	}

    if(m_savePath == "") 
	{
      if (defaultPath.isEmpty())
        defaultPath = QDir::homePath();
    } 
	
	else  defaultPath = m_savePath;

    QString defaultFileName = defaultPath + '/' + defaultFile.baseName() + "." + ext;

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
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

  // convert a first and a second letters of a GAFF atom type character array,
  // two or four letters of a UFF atom type character array,
  // and four letters of a Ghemical atom type character array 
  // into integer-type ascii code
  int PsfPotDialog::chrs2int( const char* chrs )
  {
	int retval = 0;

    char fstchr[32], scdchr[32], trdchr[32], fthchr[32];
	int fstval, scdval, trdval, fthval; 

	if( strlen(chrs) <= 2 && chrs[1] != '_' ) // GAFF, MMFF? and Partially UFF Atom Types
	{
	  sprintf( fstchr, "%i", chrs[0] );
      sprintf( scdchr, "%i", chrs[1] ); // chrs[1]=='\0' case may occur.

	  fstval = atoi(fstchr);
	  scdval = atoi(scdchr);

	  retval = fstval*100 + scdval;
	}

	else if( strlen(chrs)  <= 3 && chrs[1] == '_' ) // UFF Atom Type (case 1)
	{
	   sprintf( fstchr, "%i", chrs[0] );
       sprintf( scdchr, "%i", chrs[2] ); // chrs[2]=='\0' case may occur.

	  fstval = atoi(fstchr);
	  scdval = atoi(scdchr);

	  retval = fstval*100 + scdval;
	}

	else if( ( strlen(chrs) == 5 ) // UFF Atom Type (case 2)
		     || ( strlen(chrs) == 4 ) ) // Ghemical Atom Type
	{
	   sprintf( fstchr, "%i", chrs[0] );
       sprintf( scdchr, "%i", chrs[1] );
	   sprintf( trdchr, "%i", chrs[3] );
	   sprintf( fthchr, "%i", chrs[4] );

	  fstval = atoi(fstchr);
	  scdval = atoi(scdchr);
	  trdval = atoi(trdchr);
	  fthval = atoi(fthchr);

	  retval = fstval*1000 + scdval*100 + trdval*10 + fthval ;
	}

	return retval;
  }


  bool PsfPotDialog::bondTypeComp( const bondKind& left, const bondKind& right )
  {
	int left1, left2, right1, right2;

	left1 = chrs2int(left.aT1);
	left2 = chrs2int(left.aT2);
	right1 = chrs2int(right.aT1);
	right2 = chrs2int(right.aT2);

	if( left1 < right1 ) return true;
	else if ( left1 == right1 )
	{
		 if( left2 < right2 ) return true;
		 else return false;
	}

	else return false;	
  }

  bool PsfPotDialog::bondTypeEqual(const bondKind& left, const bondKind& right )
  {
	if ( (strcmp( left.aT1, right.aT1 ) == 0) && (strcmp( left.aT2, right.aT2 ) == 0 ) )
		return true;

	else return false;
  }

  bool PsfPotDialog::angleTypeComp( const angleKind& left, const angleKind& right )
  {
	int left1, left2, left3, right1, right2, right3;

	left1 = chrs2int(left.aT1);
	left2 = chrs2int(left.aT2);
	left3 = chrs2int(left.aT3);
	right1 = chrs2int(right.aT1);
	right2 = chrs2int(right.aT2);
	right3 = chrs2int(right.aT3);

	if( left1 < right1 ) return true;
	else if ( left1 == right1 )
	{
		 if( left2 < right2 ) return true;

		 else if( left2 == right2 )
		 {
			 if( left3 < right3 ) return true;
			 else return false;
		 }

		 else return false;
	}

	else return false;	
  }

  bool PsfPotDialog::angleTypeEqual( const angleKind& left, const angleKind& right )
  {
	if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT2 ) == 0)
		&& (strcmp( left.aT3, right.aT3 ) == 0) )
		return true;

	// swap first atom and third atom
	else if( ( strcmp( left.aT1, right.aT3 ) == 0 ) && (strcmp( left.aT2, right.aT2 ) == 0)
		&& (strcmp( left.aT3, right.aT1 ) == 0) )
		return true;

	else return false;
  }


  bool PsfPotDialog::torTypeComp( const torKind& left, const torKind& right )
  {
	int left1, left2, left3, left4;
	int right1, right2, right3, right4;

	left1 = chrs2int(left.aT1);
	left2 = chrs2int(left.aT2);
	left3 = chrs2int(left.aT3);
	left4 = chrs2int(left.aT4);
	right1 = chrs2int(right.aT1);
	right2 = chrs2int(right.aT2);
	right3 = chrs2int(right.aT3);
	right4 = chrs2int(right.aT4);

    if( left1 < right1 ) return true;

	/*
	else if ( left1 == right1 )
	{
		 if( left2 < right2 ) return true;

		 else if( left2 == right2 )
		 {
			 if( left3 < right3 ) return true;

			 else if( left3 == right3 )
			 {
				 if( left4 < right4 )
					 return true;

				 else return false;
			 }
			 else return false;
		 }

		 else return false;
	}
    */

    else if ( left1 == right1 )
	{
		 if( left2 < right2 ) return true;

		 else if( left2 == right2 )
		 {
			 if( left3 < right3 ) return true;

			 else if( left3 == right3 )
			 {
				 if( left4 < right4 )
					 return true;

				 else if( left4 == right4 )
				 {
					 if( left.gamma < right.gamma ) return true;

					 else if( left.gamma == right.gamma ) // float == float
					 {
						 if( left.n < right.n ) return true;
						 else return false;
					 }

					 else return false; 
				 }

				 else return false;
	
			 }
			 else return false;
		 }

		 else return false;
	}

	else return false;
  }


  bool PsfPotDialog::torTypeEqual( const torKind& left, const torKind& right )
  {
	/*
	if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT2 ) == 0)
		&& (strcmp( left.aT3, right.aT3 ) == 0) && (strcmp( left.aT4, right.aT4 ) == 0) )
		return true;

	// reverse order (1-2-3-4 to 4-3-2-1)
	else if( ( strcmp( left.aT1, right.aT4 ) == 0 ) && (strcmp( left.aT2, right.aT3 ) == 0)
		&& (strcmp( left.aT3, right.aT2 ) == 0) && (strcmp( left.aT4, right.aT1 ) == 0) )
		return true;

	else return false;
	*/

	if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT2 ) == 0)
		&& (strcmp( left.aT3, right.aT3 ) == 0) && (strcmp( left.aT4, right.aT4 ) == 0) 
		&& fabs(left.gamma - right.gamma) < EPS && fabs(left.n - right.n ) < EPS )
		return true;

	// reverse order (1-2-3-4 to 4-3-2-1)
	else if( ( strcmp( left.aT1, right.aT4 ) == 0 ) && (strcmp( left.aT2, right.aT3 ) == 0)
		&& (strcmp( left.aT3, right.aT2 ) == 0) && (strcmp( left.aT4, right.aT1 ) == 0) 
	    && fabs(left.gamma - right.gamma) < EPS && fabs(left.n - right.n ) < EPS )
		return true;

	else return false;

  }

  bool PsfPotDialog::oopTypeComp( const oopKind& left, const oopKind& right )
  {
	int left1, left2, left3, left4;
	int right1, right2, right3, right4;

	left1 = chrs2int(left.aT1);
	left2 = chrs2int(left.aT2);
	left3 = chrs2int(left.aT3);
	left4 = chrs2int(left.aT4);
	right1 = chrs2int(right.aT1);
	right2 = chrs2int(right.aT2);
	right3 = chrs2int(right.aT3);
	right4 = chrs2int(right.aT4);

    if( left1 < right1 ) return true;

	else if ( left1 == right1 )
	{
		 if( left2 < right2 ) return true;

		 else if( left2 == right2 )
		 {
			 if( left3 < right3 ) return true;

			 else if( left3 == right3 )
			 {
				 if( left4 < right4 )
					 return true;

				 else return false;
			 }
			 else return false;
		 }

		 else return false;
	}

	else return false;
  }


  bool PsfPotDialog::oopTypeEqual( const oopKind& left, const oopKind& right )
  {
	if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT2 ) == 0)
		&& (strcmp( left.aT3, right.aT3 ) == 0) && (strcmp( left.aT4, right.aT4 ) == 0) )
		return true;

	// permutate 2-3-4 to 2-4-3 (lexicologial order)
	else if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT2 ) == 0)
		&& (strcmp( left.aT3, right.aT4 ) == 0) && (strcmp( left.aT4, right.aT3 ) == 0) )
		return true;

	// permutate 2-3-4 to 3-2-4
    else if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT3 ) == 0)
		&& (strcmp( left.aT3, right.aT2 ) == 0) && (strcmp( left.aT4, right.aT4 ) == 0) )
		return true;

	// permutate 2-3-4 to 3-4-2
	else if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT3 ) == 0)
		&& (strcmp( left.aT3, right.aT4 ) == 0) && (strcmp( left.aT4, right.aT2 ) == 0) )
		return true;

	// permutate 2-3-4 t0 4-2-3
    else if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT4 ) == 0)
		&& (strcmp( left.aT3, right.aT2 ) == 0) && (strcmp( left.aT4, right.aT3 ) == 0) )
		return true;

	// permutate 2-3-4 to 4-3-2
    else if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT4 ) == 0)
		&& (strcmp( left.aT3, right.aT3 ) == 0) && (strcmp( left.aT4, right.aT2 ) == 0) )
		return true;

	else return false;
  }


  bool PsfPotDialog::nbTypeComp( const nbKind& left, const nbKind& right )
  {
	if( chrs2int(left.aT) < chrs2int(right.aT) ) return true;

	else return false;
  }

  bool PsfPotDialog::nbTypeEqual( const nbKind& left, const nbKind& right )
  {
    if( chrs2int(left.aT) == chrs2int(right.aT) ) return true;

	 else return false;
  }


 

  int PsfPotDialog::connAtoms( OBAtom* a )
  {
	  if( a == NULL ) return -1;

	  int ret = 0;

      FOR_NBORS_OF_ATOM( nbr, a ) ret++;

	  return ret;

  }


  	QString PsfPotDialog::psfGaff()
	{
	  QString buffer;
      QTextStream ts(&buffer);

	  QString molBaseName;
	  OBMol mol;


	  if( m_molecule != NULL )
	  {
		  	QFileInfo molNameInfo( m_molecule->fileName() );
	        molBaseName =  molNameInfo.baseName();

			mol = m_molecule->OBMol();

			// if( &mol == NULL ) return "";

			if( mol.Empty() ) 
	        { 
		       // qDebug() << "moleclar data is empty!";
		       return tr("! moleclar data is empty!");
	        }

	       // OBForceField* pFF = OBForceField::FindForceField("GAFF");
		   OBForceField* pFF = OBForceField::FindForceField("GAFF");
           if (!pFF)
	       {
               // qDebug() << "Could not find GAFF!"; 
               return tr("! Could not find GAFF!");
	       }

            if (!pFF->Setup(mol))
	       {
			   QString debugInfo;
	           QTextStream dits(&debugInfo);// debugInfo text stream;

               // qDebug() <<" Could not setup GAFF for the molecule." ;
			   dits << tr("! Could not setup GAFF for the molecule.\n");
			   dits << tr("! It may go well if you add hydrogens to the moleule.\n");
               return debugInfo;
            }

	         if(!pFF->GetAtomTypes(mol))
           {
	           // qDebug() << "Could not find atom types for the molecule." ;
	          return tr("! Could not find atom types for the molecule.");
            }

            // For multiple disconnected fragments case.
			vector<OBMol> mols = mol.Separate(1); 
			int molid;

			for ( molid = 0; molid < mols.size() ; molid++ )
			{
				FOR_ATOMS_OF_MOL( atom, mols[molid] )
				{
					char segIdx[ 1024 ];
					sprintf( segIdx, "%d", molid + 1 );
					OBPairData *segIdxLabel = new OBPairData;
                    segIdxLabel->SetAttribute("SegmentIndex");
                    segIdxLabel->SetValue(segIdx);
                      // segIdxLabel->SetOrigin(userInput); // set by user, not by Open Babel

					mol.GetAtomById( atom->GetId() )->SetData(segIdxLabel);
				}
			}


			// Set Residue Names
			/*
		    FOR_ATOMS_OF_MOL( atom, mol)
			{
		       QList<Residue*> residues = m_molecule->residues();

			   foreach (Residue *res, residues) 
		       {
		          foreach( unsigned long uid, res->atoms()) // uid is unique id of atom.  
				  {
				      if( mol.GetAtomById(uid) != NULL ) // Argument of OBMol::GetAtomById(arg) is also unique id.
				      {
				         OBResidue *tempRes = mol.GetAtomById(uid)->GetResidue();
			             tempRes->SetName( res->name().toStdString() );
			             mol.GetAtomById(uid)->SetResidue(tempRes);
				       }
			       }
		        }

			}
			*/

			ts << "PSF\n\n";
            ts << "      3 !NTITLE\n";
			ts << "      " << left << molBaseName << ".\n";
            ts << tr("      Avogadro generated Protein Structure File (PSF)\n");
			ts << tr("      by using General AMBER Force Fields (GAFF).\n\n");

			ts << qSetFieldWidth(5) << right << mol.NumAtoms() << " !NATOMS" << qSetFieldWidth(0) << endl;


			  FOR_ATOMS_OF_MOL( atom, mol )
		     {
				 int atomid = atom->GetIdx();

			      OBPairData *type = (OBPairData*) atom->GetData("FFAtomType");
                  OBPairData *segidx;
				  if( atom->HasData("SegmentIndex"))
				      segidx =  (OBPairData*) atom->GetData("SegmentIndex");

				  if( atom->HasData("SegmentIndex"))
				  { ts << qSetFieldWidth(7) << right << atomid << qSetFieldWidth(0)
				 	 << " SEG" << left << segidx->GetValue().c_str();
				  }

				  else
				  {
					  ts << qSetFieldWidth(7) << right << atomid << qSetFieldWidth(0)
				 	      << " SEG" ;
				  }

				  ts << qSetFieldWidth(7) << right 
					  << atom->GetResidue()->GetIdx() + 1;
				  ts << qSetFieldWidth(5) << right << atom->GetResidue()->GetName().c_str() 
					  << qSetFieldWidth(0) << " ";
				  ts << left << etab.GetSymbol(atom->GetAtomicNum()) << "  ";
				  ts << left << type->GetValue().c_str() << "    ";
				  ts << qSetRealNumberPrecision(5) << fixed << right << atom->GetPartialCharge() <<  "  ";
				  ts << qSetRealNumberPrecision(5) << fixed << right << atom->GetExactMass() << "      0\n";

			  }	   

	      ts << endl;
          ts << qSetFieldWidth(5) << right << mol.NumBonds() << " !NBONDS" << qSetFieldWidth(0) << endl;
		  ts << "     ";

		 // bond pairs
         int i = 0;
		 
		 FOR_BONDS_OF_MOL( bond, mol )
          {
				ts << "   " 
				    << bond->GetBeginAtom()->GetIdx()
				    << "   " 
				   << bond->GetEndAtom()->GetIdx();
	               i++;

	              if( i >= 4 ) // four pairs of atoms per line
	              {
			             ts << endl << "     ";
	                     i=0;
	              }
	              else continue;
           }
          // end of bonds

		  ts << endl << endl;

		  // number of angles
          int nAngle = 0;
          FOR_ANGLES_OF_MOL(angle, mol){ nAngle++; }

          ts << qSetFieldWidth(5) << right << nAngle << " !NTHETA: angles" << qSetFieldWidth(0) << endl;
		  ts << "     ";

		  // angle triples
          OBAtom *a, *b, *c;
          i = 0;

		   	    FOR_ANGLES_OF_MOL( angle, mol )
		        {					
		             a = mol.GetAtom((*angle)[1] + 1);
	                 b = mol.GetAtom((*angle)[0] + 1); // vertex
	                 c = mol.GetAtom((*angle)[2] + 1);

	                ts << "   " 
			        	<< a->GetIdx()
			    	    << "   " 
				        << b->GetIdx()
				         << "   " 
				       << c->GetIdx();
                     i++;

	                if( i >= 3 ) // three triples of atoms per line
	                { 
		                ts << endl << "     ";
		                 i = 0;
	                }
	                else continue;

		         }
          ts << endl;
         // end of angles

		  ts << endl; // an empty line

        // number of torsion (dihedral) quadreples
        int nTorsion = 0;
        FOR_TORSIONS_OF_MOL(torsion, mol){ nTorsion++ ;}
		ts << qSetFieldWidth(5) << right << nTorsion << " !NPHI: dihedrals" << qSetFieldWidth(0) << endl;
	    ts << "     ";

        // torsion (dihedral) quadreples
        OBAtom* d; // OBAtom *a, *b and *c are already decleard above.
        i = 0;

	     FOR_TORSIONS_OF_MOL(torsion, mol )
         {
	         a = mol.GetAtom((*torsion)[0] + 1);
	         b = mol.GetAtom((*torsion)[1] + 1);
	         c = mol.GetAtom((*torsion)[2] + 1);
	         d = mol.GetAtom((*torsion)[3] + 1);

			 ts << "   " 
				  << a->GetIdx()
				  << "   " 
				  << b->GetIdx()
				  << "   "
				  << c->GetIdx()
				  << "   " <<
				     d->GetIdx();
	           i++;

	          if( i >= 2 ) // two quaruples of atoms per line
	          { 
		        ts << endl
					<< "     ";
		        i = 0;
	           }
	        else continue;
         }
         // end of dihedrals

          ts << endl
		  	 << "     "; // an empty line
 
         // number of out-of-plane ( improper ) dihedral quadreples
          int nOOP = 0;
          FOR_ATOMS_OF_MOL (atom, mol) 
         { 
			 if(  connAtoms(&*atom) == 3 ) {nOOP++;}
         }
         ts << nOOP << " !NIMPHI: impropers" << qSetFieldWidth(0) << endl;
		 ts << "     ";

        // out-of-plane (improper) dihedral quadreples
        i = 0;

		FOR_ATOMS_OF_MOL ( atom, mol ) 
        {
			 if(  connAtoms(&*atom) == 3 )
	         {
				  ts << "   " << atom->GetIdx();

                   FOR_NBORS_OF_ATOM( nbr, &*atom )
				   { ts << "   " << nbr->GetIdx(); }
		          i++;

		          if( i >= 2 ) // two quaruples of atoms per line
		    	  {
		             ts << endl << "     ";
				     i = 0;
		          }

			  }
	           else continue;
		 }
       // end of out-of-plane

        ts << endl; // an empty line

        ts << "   0 !NDON\n\n";
        ts << "   0 !NACC\n\n";
        ts << "   0 !NNB\n\n";
        ts << "   0 !NGRP\n\n";
	 }

	  else
	  {
		  ts << tr("! No molecule was found.\n");
		  ts << tr("! These psf and pot preview panes are NOT automatically updated !\n");
		  ts << tr("! Please click the update button to update !\n");
	  }

	  return buffer;
	}
 // psfGaff() ends.


	QString PsfPotDialog::potGaff()
	{
	  QString buffer;
	  QTextStream ts(&buffer);

	  QString molBaseName;
	  OBMol mol;

	  if( m_molecule != NULL )
	  {
		   QFileInfo molNameInfo( m_molecule->fileName() );
	        molBaseName =  molNameInfo.baseName();

		    mol = m_molecule->OBMol();

			if( mol.Empty() ) 
	        { 
		       // qDebug() << "moleclar data is empty!";
		       return tr("! moleclar data is empty!");
	        }

	       OBForceField* pFF = OBForceField::FindForceField("GAFF");

           if (!pFF)
	       {
               // qDebug() << "Could not find GAFF!"; 
               return tr("! Could not find GAFF!");
	       }

           if (!pFF->Setup(mol))
	       {
			   QString debugInfo;
	           QTextStream dits(&debugInfo); // debugInfo text stream;

               // qDebug() <<" Could not setup GAFF for the molecule." ;
			   dits << tr("! Could not setup GAFF for the molecule.\n");
			   dits << tr("! It may go well if you add hydrogens to the moleule.\n");
               return debugInfo;
            }

	       if(!pFF->GetAtomTypes(mol))
           {
	           // qDebug() << "Could not find atom types for the molecule." ;
	          return tr("! Could not find atom types for the molecule.");
           }

		   ffpGaff* pP = static_cast<ffpGaff*>(pFF);
		   if (pP->ffpGaff::SetupPointers() == false )
		   { return tr("Failed to get GAFF Parameters");}

		   sort( bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeComp );
	       vector<bondKind>::iterator bondEnd = unique(bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeEqual );
	       bkmol.erase( bondEnd, bkmol.end() );

		   sort( akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeComp );
	       vector<angleKind>::iterator angleEnd = unique(akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeEqual );
	       akmol.erase( angleEnd, akmol.end() );

		   sort( tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeComp );
	       vector<torKind>::iterator torEnd = unique(tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeEqual );
	       tkmol.erase( torEnd, tkmol.end() );

		   sort( okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeComp );
	       vector<oopKind>::iterator oopEnd = unique(okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeEqual );
	       okmol.erase( oopEnd, okmol.end() );

		   sort( nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeComp );
	       vector<nbKind>::iterator nbEnd = unique(nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeEqual );
	       nkmol.erase( nbEnd, nkmol.end() );


		  ts << tr("*>>>>>>>   General AMBER FF (GAFF) written in CHARMM FF style  <<<<<<<\n");
		  ts << tr("*>>>>>>>   for ") << molBaseName;
		  ts << tr(" which generated by Avogadro.                         <<<<<<<\n\n");
		  ts << "BONDS\n";
		  ts << "!\n";
          ts << "!V(bond) = Kb(b - b0)**2\n";
          ts << "!\n";
          ts << "!Kb: kcal/mole/A**2\n";
          ts << "!b0: A\n";
          ts << "!\n";
          ts << "!atom type Kb          b0\n";
          ts << "!\n";

		 
		 for( int i = 0; i < bkmol.size(); i++ )
	     {

		  ts << qSetRealNumberPrecision(5) << fixed << right
			 << bkmol[i].aT1 << " "
             <<  bkmol[i].aT2 << "    "
			 << qSetRealNumberPrecision(9)
			 <<  bkmol[i].Kb << " "
			 <<  bkmol[i].req << endl;
   	     }
		 

		  ts << endl; // insert empty line

		 ts << "ANGLES\n";
		 ts << "!\n";
		 ts << "!V(angle) = Ktheta(Theta - Theta0)**2\n";
		 ts << "!\n";
		 ts << "!V(Urey-Bradley) = Kub(S - S0)**2\n";
		 ts << "!\n";
		 ts << "!Ktheta: kcal/mole/rad**2\n";
		 ts << "!Theta0: degrees\n";
		 ts << "!Kub: kcal/mole/A**2 (Urey-Bradley)\n";
		 ts << "!S0: A\n";
		 ts << "!\n";
		 ts << tr("!atom types     Ktheta    Theta0   Kub     S0\n");
		 ts << "!\n";

		 for( int i = 0; i < akmol.size(); i++ )
	     {
		    ts  << qSetRealNumberPrecision(5) << fixed << right
				<<  akmol[i].aT1 << " "
               <<  akmol[i].aT2 << " "
			   <<  akmol[i].aT3 << "    "
			   <<  qSetRealNumberPrecision(9)
			   <<  akmol[i].Kth << " "
			   <<  akmol[i].theq << endl;
	     
		 }


		ts << endl; // an empty line

	    ts << "DIHEDRALS" << endl;
        ts << "!\n";
        ts << "!V(dihedral) = Kchi(1 + cos(n(chi) - delta))\n";
        ts << "!\n";
        ts << "!Kchi: kcal/mole\n";
        ts << "!n: multiplicity\n";
        ts << "!delta: degrees\n";
        ts << "!\n";
        ts << "!atom types             Kchi    n   delta\n";
		ts << "!\n";


	   for( int i = 0; i < tkmol.size(); i++ )
	   {
		      ts  << qSetRealNumberPrecision(5) << fixed << right
			      <<  tkmol[i].aT1 << " "
                  <<  tkmol[i].aT2 << " "
			      <<  tkmol[i].aT3 << " "
                  <<  tkmol[i].aT4 << "    "
			      <<  tkmol[i].vn2 << " "
			      <<  qSetRealNumberPrecision(0) 
			      <<  tkmol[i].n << "  "
			      <<  qSetRealNumberPrecision(5) << fixed << right
			      <<  tkmol[i].gamma << endl;
   	    }


	    ts << endl; // an empty line

	    ts << "IMPROPER" << endl;
        ts << "!V(improper) = Kpsi(psi - psi0)**2\n";
        ts << "!\n";
        ts << "!Kpsi: kcal/mole/rad**2\n";
        ts << "!psi0: degrees\n";
        ts << "!note that the second column of numbers (0) is ignored\n";
        ts << "!\n";


	   for( int i = 0; i < okmol.size(); i++ )
	  {
		    ts << qSetRealNumberPrecision(5) << fixed << right
		    	<<  okmol[i].aT1 << " "
               <<  okmol[i].aT2 << " "
			   <<  okmol[i].aT3 << " "
               <<  okmol[i].aT4 << "    "
			   <<  okmol[i].vn2 << "   "
			   << "0    "
			   <<  okmol[i].gamma << endl;

   	  }


	ts << endl; // an empty line

	ts << "NONBONDED" << endl;
	ts << "!\n";
	ts << "!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]\n";
	ts << "!\n";
	ts << "!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)\n";
	ts << "!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j\n";
	ts << "!\n";
	ts << "!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4\n";
	ts << "!\n";


	for( int i = 0; i< nkmol.size() ; i++ )
	{

		ts  << qSetRealNumberPrecision(5) << fixed << right
			<< nkmol[i].aT << "     "
			<< "0.00000" << "   "
			<< nkmol[i].eps <<  "   "
			<< nkmol[i].R
			<< endl;

	}


	ts << endl; // an empty line

	ts << "END" << endl << endl; 

	  }// end( if(m_molecule != NULL)

	  else
	  {
		  ts << tr("! No molecule was found.\n");
		  ts << tr("! These psf and pot preview panes are NOT automatically updated !\n");
		  ts << tr("! Please click the update button to update !\n");
	  }

	  return buffer;

	} 


	QString PsfPotDialog::psfUff()
	{
	  QString buffer;
	  QTextStream ts(&buffer);

	  // ts << tr("UFF is not yet implemented for psf.") << endl;

	  QString molBaseName;
	  OBMol mol;


	  if( m_molecule != NULL )
	  {
		  	QFileInfo molNameInfo( m_molecule->fileName() );
	        molBaseName =  molNameInfo.baseName();

			mol = m_molecule->OBMol();

			if( mol.Empty() ) 
	        { 
		       return tr("! moleclar data is empty!");
	        }

		   OBForceField* pFF = OBForceField::FindForceField("UFF");
           if (!pFF)
	       {
               return tr("! Could not find UFF!");
	       }

            if (!pFF->Setup(mol))
	       {
			   QString debugInfo;
	           QTextStream dits(&debugInfo);// debugInfo text stream;

			   dits << tr("! Could not setup UFF for the molecule.\n");
			   dits << tr("! It may go well if you add hydrogens to the moleule.\n");
               return debugInfo;
            }

	        if(!pFF->GetAtomTypes(mol))
           {
	          return tr("! Could not find atom types for the molecule.");
           }

			// For multiple disconnected fragments case.
			vector<OBMol> mols = mol.Separate(1); 
			int molid;

			for ( molid = 0; molid < mols.size() ; molid++ )
			{
				FOR_ATOMS_OF_MOL( atom, mols[molid] )
				{
					char segIdx[ 1024 ];
					sprintf( segIdx, "%d", molid + 1 );
					OBPairData *segIdxLabel = new OBPairData;
                    segIdxLabel->SetAttribute("SegmentIndex");
                    segIdxLabel->SetValue(segIdx);
                      // segIdxLabel->SetOrigin(userInput); // set by user, not by Open Babel

					mol.GetAtomById( atom->GetId() )->SetData(segIdxLabel);
				}
			}

			ts << "PSF\n\n";
            ts << "      3 !NTITLE\n";
			ts << "      " << left << molBaseName << ".\n";
            ts << tr("      Avogadro generated Protein Structure File (PSF)\n");
			ts << tr("      by using Universal Force Fields (UFF).\n\n");

			ts << qSetFieldWidth(5) << right << mol.NumAtoms() << " !NATOMS" << qSetFieldWidth(0) << endl;


			  FOR_ATOMS_OF_MOL( atom, mol )
		     {
				 int atomid = atom->GetIdx();

			      OBPairData *type = (OBPairData*) atom->GetData("FFAtomType");
                  OBPairData *segidx;
				  if( atom->HasData("SegmentIndex"))
				      segidx =  (OBPairData*) atom->GetData("SegmentIndex");

				  if( atom->HasData("SegmentIndex"))
				  { ts << qSetFieldWidth(7) << right << atomid << qSetFieldWidth(0)
				 	 << " SEG" << left << segidx->GetValue().c_str();
				  }

				  else
				  {
					  ts << qSetFieldWidth(7) << right << atomid << qSetFieldWidth(0)
				 	      << " SEG" ;
				  }

				  ts << qSetFieldWidth(7) << right 
					  << atom->GetResidue()->GetIdx() + 1;
				  ts << qSetFieldWidth(5) << right << atom->GetResidue()->GetName().c_str() 
					  << qSetFieldWidth(0) << " ";
				  ts << left << etab.GetSymbol(atom->GetAtomicNum()) << "  ";
				  ts << left << type->GetValue().c_str() << "    ";
				  ts << qSetRealNumberPrecision(5) << fixed << right << atom->GetPartialCharge() <<  "  ";
				  ts << qSetRealNumberPrecision(5) << fixed << right << atom->GetExactMass() << "      0\n";

			  }	   


			  	      ts << endl;
          ts << qSetFieldWidth(5) << right << mol.NumBonds() << " !NBONDS" << qSetFieldWidth(0) << endl;
		  ts << "     ";

		 // bond pairs
         int i = 0;
		 
		 FOR_BONDS_OF_MOL( bond, mol )
          {
				ts << "   " 
				    << bond->GetBeginAtom()->GetIdx()
				    << "   " 
				   << bond->GetEndAtom()->GetIdx();
	               i++;

	              if( i >= 4 ) // four pairs of atoms per line
	              {
			             ts << endl << "     ";
	                     i=0;
	              }
	              else continue;
           }
          // end of bonds

		  ts << endl << endl;

		  // number of angles
          int nAngle = 0;
          FOR_ANGLES_OF_MOL(angle, mol){ nAngle++; }

          ts << qSetFieldWidth(5) << right << nAngle << " !NTHETA: angles" << qSetFieldWidth(0) << endl;
		  ts << "     ";

		  // angle triples
          OBAtom *a, *b, *c;
          i = 0;

		   	    FOR_ANGLES_OF_MOL( angle, mol )
		        {					
		             a = mol.GetAtom((*angle)[1] + 1);
	                 b = mol.GetAtom((*angle)[0] + 1); // vertex
	                 c = mol.GetAtom((*angle)[2] + 1);

	                ts << "   " 
			        	<< a->GetIdx()
			    	    << "   " 
				        << b->GetIdx()
				         << "   " 
				       << c->GetIdx();
                     i++;

	                if( i >= 3 ) // three triples of atoms per line
	                { 
		                ts << endl << "     ";
		                 i = 0;
	                }
	                else continue;

		         }
          ts << endl;
         // end of angles

		  ts << endl; // an empty line

        // number of torsion (dihedral) quadreples
        int nTorsion = 0;
        FOR_TORSIONS_OF_MOL(torsion, mol){ nTorsion++ ;}
		ts << qSetFieldWidth(5) << right << nTorsion << " !NPHI: dihedrals" << qSetFieldWidth(0) << endl;
	    ts << "     ";

        // torsion (dihedral) quadreples
        OBAtom* d; // OBAtom *a, *b and *c are already decleard above.
        i = 0;

	     FOR_TORSIONS_OF_MOL(torsion, mol )
         {
	         a = mol.GetAtom((*torsion)[0] + 1);
	         b = mol.GetAtom((*torsion)[1] + 1);
	         c = mol.GetAtom((*torsion)[2] + 1);
	         d = mol.GetAtom((*torsion)[3] + 1);

			 ts << "   " 
				  << a->GetIdx()
				  << "   " 
				  << b->GetIdx()
				  << "   "
				  << c->GetIdx()
				  << "   " <<
				     d->GetIdx();
	           i++;

	          if( i >= 2 ) // two quaruples of atoms per line
	          { 
		        ts << endl
					<< "     ";
		        i = 0;
	           }
	        else continue;
         }
         // end of dihedrals

          ts << endl
		  	 << "     "; // an empty line
 
         // number of out-of-plane ( improper ) dihedral quadreples
          int nOOP = 0;
          FOR_ATOMS_OF_MOL (atom, mol) 
         { 
			 if(  connAtoms(&*atom) == 3 ) {nOOP++;}
         }
         ts << nOOP << " !NIMPHI: impropers" << qSetFieldWidth(0) << endl;
		 ts << "     ";

        // out-of-plane (improper) dihedral quadreples
        i = 0;

		FOR_ATOMS_OF_MOL ( atom, mol ) 
        {
			 if(  connAtoms(&*atom) == 3 )
	         {
				  ts << "   " << atom->GetIdx();

                   FOR_NBORS_OF_ATOM( nbr, &*atom )
				   { ts << "   " << nbr->GetIdx(); }
		          i++;

		          if( i >= 2 ) // two quaruples of atoms per line
		    	  {
		             ts << endl << "     ";
				     i = 0;
		          }

			  }
	           else continue;
		 }
       // end of out-of-plane

        ts << endl; // an empty line

        ts << "   0 !NDON\n\n";
        ts << "   0 !NACC\n\n";
        ts << "   0 !NNB\n\n";
        ts << "   0 !NGRP\n\n";
	 }

	  else
	  {
		  ts << tr("! No molecule was found.\n");
		  ts << tr("! These psf and pot preview panes are NOT automatically updated !\n");
		  ts << tr("! Please click the update button to update !\n");
	  }

	  return buffer;

	}// psfUff() ends.


    QString PsfPotDialog::potUff()
	{
	  QString buffer;
	  QTextStream ts(&buffer);

	  // ts << tr("UFF is not yet implemented for pot.") << endl;

	  QString molBaseName;
	  OBMol mol;

	  if( m_molecule != NULL )
	  {
		   QFileInfo molNameInfo( m_molecule->fileName() );
	        molBaseName =  molNameInfo.baseName();

		    mol = m_molecule->OBMol();

			if( mol.Empty() ) 
	        { 
		       return tr("! moleclar data is empty!");
	        }

	       OBForceField* pFF = OBForceField::FindForceField("UFF");

           if (!pFF)
	       {
               return tr("! Could not find UFF!");
	       }

           if (!pFF->Setup(mol))
	       {
			   QString debugInfo;
	           QTextStream dits(&debugInfo); // debugInfo text stream;

			   dits << tr("! Could not setup UFF for the molecule.\n");
			   dits << tr("! It may go well if you add hydrogens to the moleule.\n");
               return debugInfo;
            }

	       if(!pFF->GetAtomTypes(mol))
           {
	          return tr("! Could not find atom types for the molecule.");
           }

		   ffpUff* pP = static_cast<ffpUff*>(pFF);
		   if (pP->ffpUff::SetupPointers() == false )
		   { return tr("Failed to get UFF Parameters");}

		   sort( bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeComp );
	       vector<bondKind>::iterator bondEnd = unique(bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeEqual );
	       bkmol.erase( bondEnd, bkmol.end() );

		   sort( akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeComp );
	       vector<angleKind>::iterator angleEnd = unique(akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeEqual );
	       akmol.erase( angleEnd, akmol.end() );

		   sort( tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeComp );
	       vector<torKind>::iterator torEnd = unique(tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeEqual );
	       tkmol.erase( torEnd, tkmol.end() );

		   sort( okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeComp );
	       vector<oopKind>::iterator oopEnd = unique(okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeEqual );
	       okmol.erase( oopEnd, okmol.end() );

		   sort( nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeComp );
	       vector<nbKind>::iterator nbEnd = unique(nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeEqual );
	       nkmol.erase( nbEnd, nkmol.end() );

		  ts << tr("*>>>>>>>   Universal Force Field (UFF) written in CHARMM FF style  <<<<<<<\n");
		  ts << tr("*>>>>>>>   for ") << molBaseName;
		  ts << tr(" which generated by Avogadro.                         <<<<<<<\n\n");
		  ts << "BONDS\n";
		  ts << "!\n";
          ts << "!V(bond) = Kb(b - b0)**2\n";
          ts << "!\n";
          ts << "!Kb: kcal/mole/A**2\n";
          ts << "!b0: A\n";
          ts << "!\n";
          ts << "!atom type Kb          b0\n";
          ts << "!\n";
		 
		 for( int i = 0; i < bkmol.size(); i++ )
	     {

		  ts << qSetRealNumberPrecision(5) << fixed << right
			 << bkmol[i].aT1 << " "
             <<  bkmol[i].aT2 << "    "
			 << qSetRealNumberPrecision(9)
			 <<  bkmol[i].Kb << " "
			 <<  bkmol[i].req << endl;
   	     }
		 
		  ts << endl; // insert empty line

		 ts << "ANGLES\n";
		 ts << "!\n";
		 ts << "!V(angle) = Ktheta(Theta - Theta0)**2\n";
		 ts << "!\n";
		 ts << "!V(Urey-Bradley) = Kub(S - S0)**2\n";
		 ts << "!\n";
		 ts << "!Ktheta: kcal/mole/rad**2\n";
		 ts << "!Theta0: degrees\n";
		 ts << "!Kub: kcal/mole/A**2 (Urey-Bradley)\n";
		 ts << "!S0: A\n";
		 ts << "!\n";
		 ts << tr("!atom types     Ktheta    Theta0   Kub     S0\n");
		 ts << "!\n";

		 for( int i = 0; i < akmol.size(); i++ )
	     {
		    ts  << qSetRealNumberPrecision(5) << fixed << right
				<<  akmol[i].aT1 << " "
               <<  akmol[i].aT2 << " "
			   <<  akmol[i].aT3 << "    "
			   <<  qSetRealNumberPrecision(9)
			   <<  akmol[i].Kth << " "
			   <<  akmol[i].theq << endl;
	     
		 }

		ts << endl; // an empty line

	    ts << "DIHEDRALS" << endl;
        ts << "!\n";
        ts << "!V(dihedral) = Kchi(1 + cos(n(chi) - delta))\n";
        ts << "!\n";
        ts << "!Kchi: kcal/mole\n";
        ts << "!n: multiplicity\n";
        ts << "!delta: degrees\n";
        ts << "!\n";
        ts << "!atom types             Kchi    n   delta\n";
		ts << "!\n";


	   for( int i = 0; i < tkmol.size(); i++ )
	   {
		      ts  << qSetRealNumberPrecision(5) << fixed << right
			      <<  tkmol[i].aT1 << " "
                  <<  tkmol[i].aT2 << " "
			      <<  tkmol[i].aT3 << " "
                  <<  tkmol[i].aT4 << "    "
			      <<  tkmol[i].vn2 << " "
			      <<  qSetRealNumberPrecision(0) 
			      <<  tkmol[i].n << "  "
			      <<  qSetRealNumberPrecision(5) << fixed << right
			      <<  tkmol[i].gamma << endl;
   	    }


	    ts << endl; // an empty line

	    ts << "IMPROPER" << endl;
        ts << "!V(improper) = Kpsi(psi - psi0)**2\n";
        ts << "!\n";
        ts << "!Kpsi: kcal/mole/rad**2\n";
        ts << "!psi0: degrees\n";
        ts << "!note that the second column of numbers (0) is ignored\n";
        ts << "!\n";


	   for( int i = 0; i < okmol.size(); i++ )
	  {
		    ts << qSetRealNumberPrecision(5) << fixed << right
		    	<<  okmol[i].aT1 << " "
               <<  okmol[i].aT2 << " "
			   <<  okmol[i].aT3 << " "
               <<  okmol[i].aT4 << "    "
			   <<  okmol[i].vn2 << "   "
			   << "0    "
			   <<  okmol[i].gamma << endl;

   	  }


	ts << endl; // an empty line

	ts << "NONBONDED" << endl;
	ts << "!\n";
	ts << "!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]\n";
	ts << "!\n";
	ts << "!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)\n";
	ts << "!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j\n";
	ts << "!\n";
	ts << "!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4\n";
	ts << "!\n";


	for( int i = 0; i< nkmol.size() ; i++ )
	{

		ts  << qSetRealNumberPrecision(5) << fixed << right
			<< nkmol[i].aT << "     "
			<< "0.00000" << "   "
			<< nkmol[i].eps <<  "   "
			<< nkmol[i].R
			<< endl;

	}


	ts << endl; // an empty line

	ts << "END" << endl << endl; 

	  }// end( if(m_molecule != NULL)

	  else
	  {
		  ts << tr("! No molecule was found.\n");
		  ts << tr("! These psf and pot preview panes are NOT automatically updated !\n");
		  ts << tr("! Please click the update button to update !\n");
	  }

	  return buffer;

   }// end potUff()

   QString PsfPotDialog::psfGhemical()
   {
	  // return "! Ghemical psf is not yet implemented !";

	  QString buffer;
	  QTextStream ts(&buffer);

	  QString molBaseName;
	  OBMol mol;

	  if( m_molecule != NULL )
	  {
		  	QFileInfo molNameInfo( m_molecule->fileName() );
	        molBaseName =  molNameInfo.baseName();

			mol = m_molecule->OBMol();

			if( mol.Empty() ) 
	        { 
		       return tr("! moleclar data is empty!");
	        }

		   OBForceField* pFF = OBForceField::FindForceField("Ghemical");
           if (!pFF)
	       {
               return tr("! Could not find Ghemical foecefield !");
	       }

            if (!pFF->Setup(mol))
	       {
			   QString debugInfo;
	           QTextStream dits(&debugInfo);// debugInfo text stream;

			   dits << tr("! Could not setup Ghemical forcefield for the molecule.\n");
			   dits << tr("! It may go well if you add hydrogens to the moleule.\n");
               return debugInfo;
            }

	        if(!pFF->GetAtomTypes(mol))
           {
	          return tr("! Could not find atom types for the molecule.");
           }

			// For multiple disconnected fragments case.
			vector<OBMol> mols = mol.Separate(1); 
			int molid;

			for ( molid = 0; molid < mols.size() ; molid++ )
			{
				FOR_ATOMS_OF_MOL( atom, mols[molid] )
				{
					char segIdx[ 1024 ];
					sprintf( segIdx, "%d", molid + 1 );
					OBPairData *segIdxLabel = new OBPairData;
                    segIdxLabel->SetAttribute("SegmentIndex");
                    segIdxLabel->SetValue(segIdx);
                      // segIdxLabel->SetOrigin(userInput); // set by user, not by Open Babel

					mol.GetAtomById( atom->GetId() )->SetData(segIdxLabel);
				}
			}

			ts << "PSF\n\n";
            ts << "      3 !NTITLE\n";
			ts << "      " << left << molBaseName << ".\n";
            ts << tr("      Avogadro generated Protein Structure File (PSF)\n");
			ts << tr("      by using Ghemical Force Fields.\n\n");

			ts << qSetFieldWidth(5) << right << mol.NumAtoms() << " !NATOMS" << qSetFieldWidth(0) << endl;


			  FOR_ATOMS_OF_MOL( atom, mol )
		     {
				 int atomid = atom->GetIdx();

			      OBPairData *type = (OBPairData*) atom->GetData("FFAtomType");
                  OBPairData *segidx;
				  if( atom->HasData("SegmentIndex"))
				      segidx =  (OBPairData*) atom->GetData("SegmentIndex");

				  if( atom->HasData("SegmentIndex"))
				  { ts << qSetFieldWidth(7) << right << atomid << qSetFieldWidth(0)
				 	 << " SEG" << left << segidx->GetValue().c_str();
				  }

				  else
				  {
					  ts << qSetFieldWidth(7) << right << atomid << qSetFieldWidth(0)
				 	      << " SEG" ;
				  }

				  ts << qSetFieldWidth(7) << right 
					  << atom->GetResidue()->GetIdx() + 1;
				  ts << qSetFieldWidth(5) << right << atom->GetResidue()->GetName().c_str() 
					  << qSetFieldWidth(0) << " ";
				  ts << left << etab.GetSymbol(atom->GetAtomicNum()) << "  ";
				  ts << left << type->GetValue().c_str() << "    ";
				  ts << qSetRealNumberPrecision(5) << fixed << right << atom->GetPartialCharge() <<  "  ";
				  ts << qSetRealNumberPrecision(5) << fixed << right << atom->GetExactMass() << "      0\n";

			  }	   


			  	      ts << endl;
          ts << qSetFieldWidth(5) << right << mol.NumBonds() << " !NBONDS" << qSetFieldWidth(0) << endl;
		  ts << "     ";

		 // bond pairs
         int i = 0;
		 
		 FOR_BONDS_OF_MOL( bond, mol )
          {
				ts << "   " 
				    << bond->GetBeginAtom()->GetIdx()
				    << "   " 
				   << bond->GetEndAtom()->GetIdx();
	               i++;

	              if( i >= 4 ) // four pairs of atoms per line
	              {
			             ts << endl << "     ";
	                     i=0;
	              }
	              else continue;
           }
          // end of bonds

		  ts << endl << endl;

		  // number of angles
          int nAngle = 0;
          FOR_ANGLES_OF_MOL(angle, mol){ nAngle++; }

          ts << qSetFieldWidth(5) << right << nAngle << " !NTHETA: angles" << qSetFieldWidth(0) << endl;
		  ts << "     ";

		  // angle triples
          OBAtom *a, *b, *c;
          i = 0;

		   	    FOR_ANGLES_OF_MOL( angle, mol )
		        {					
		             a = mol.GetAtom((*angle)[1] + 1);
	                 b = mol.GetAtom((*angle)[0] + 1); // vertex
	                 c = mol.GetAtom((*angle)[2] + 1);

	                ts << "   " 
			        	<< a->GetIdx()
			    	    << "   " 
				        << b->GetIdx()
				         << "   " 
				       << c->GetIdx();
                     i++;

	                if( i >= 3 ) // three triples of atoms per line
	                { 
		                ts << endl << "     ";
		                 i = 0;
	                }
	                else continue;

		         }
          ts << endl;
         // end of angles

		ts << endl; // an empty line

        // number of torsion (dihedral) quadreples
        int nTorsion = 0;
        FOR_TORSIONS_OF_MOL(torsion, mol){ nTorsion++ ;}
		ts << qSetFieldWidth(5) << right << nTorsion << " !NPHI: dihedrals" << qSetFieldWidth(0) << endl;
	    ts << "     ";

        // torsion (dihedral) quadreples
        OBAtom* d; // OBAtom *a, *b and *c are already decleard above.
        i = 0;

	     FOR_TORSIONS_OF_MOL(torsion, mol )
         {
	         a = mol.GetAtom((*torsion)[0] + 1);
	         b = mol.GetAtom((*torsion)[1] + 1);
	         c = mol.GetAtom((*torsion)[2] + 1);
	         d = mol.GetAtom((*torsion)[3] + 1);

			 ts << "   " 
				  << a->GetIdx()
				  << "   " 
				  << b->GetIdx()
				  << "   "
				  << c->GetIdx()
				  << "   " <<
				     d->GetIdx();
	           i++;

	          if( i >= 2 ) // two quaruples of atoms per line
	          { 
		        ts << endl
					<< "     ";
		        i = 0;
	           }
	        else continue;
         }
       // end of dihedrals

      ts << endl; // an empty line

	  // improper (out-of plane) motion energy is ignored in ghemical ff.
       ts <<  "   0 !NIMPHI: impropers" << endl;

       // end of out-of-plane

        ts << endl; // an empty line

        ts << "   0 !NDON\n\n";
        ts << "   0 !NACC\n\n";
        ts << "   0 !NNB\n\n";
        ts << "   0 !NGRP\n\n";
	 }

	  else
	  {
		  ts << tr("! No molecule was found.\n");
		  ts << tr("! These psf and pot preview panes are NOT automatically updated !\n");
		  ts << tr("! Please click the update button to update !\n");
	  }

	  return buffer;

   }// end psfGhemical()

   QString PsfPotDialog::potGhemical()
   {
	  // return "! Ghemical pot is not yet implemented !";

	  QString buffer;
	  QTextStream ts(&buffer);

	  QString molBaseName;
	  OBMol mol;

	  if( m_molecule != NULL )
	  {
		   QFileInfo molNameInfo( m_molecule->fileName() );
	        molBaseName =  molNameInfo.baseName();

		    mol = m_molecule->OBMol();

			if( mol.Empty() ) 
	        { 
		       return tr("! moleclar data is empty!");
	        }

	       OBForceField* pFF = OBForceField::FindForceField("Ghemical");

           if (!pFF)
	       {
               return tr("! Could not find Ghemical forcefield!");
	       }

           if (!pFF->Setup(mol))
	       {
			   QString debugInfo;
	           QTextStream dits(&debugInfo); // debugInfo text stream;

			   dits << tr("! Could not setup Ghemical forcefield for the molecule.\n");
			   dits << tr("! It may go well if you add hydrogens to the moleule.\n");
               return debugInfo;
            }

	       if(!pFF->GetAtomTypes(mol))
           {
	          return tr("! Could not find atom types for the molecule.");
           }

		   ffpGhemical* pP = static_cast<ffpGhemical*>(pFF);
		   if (pP->ffpGhemical::SetupPointers() == false )
		   { return tr("Failed to get Ghemical forcefield Parameters");}


		   sort( bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeComp );
	       vector<bondKind>::iterator bondEnd = unique(bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeEqual );
	       bkmol.erase( bondEnd, bkmol.end() );



		   sort( akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeComp );
	       vector<angleKind>::iterator angleEnd = unique(akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeEqual );
	       akmol.erase( angleEnd, akmol.end() );



		   sort( tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeComp );
	       vector<torKind>::iterator torEnd = unique(tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeEqual );
	       tkmol.erase( torEnd, tkmol.end() );

           
		   /*
		   sort( okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeComp );
	       vector<oopKind>::iterator oopEnd = unique(okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeEqual );
	       okmol.erase( oopEnd, okmol.end() );
		   */


		   sort( nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeComp );
	       vector<nbKind>::iterator nbEnd = unique(nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeEqual );
	       nkmol.erase( nbEnd, nkmol.end() );


		  ts << tr("*>>>>>>>   Ghemical Force Field written in CHARMM FF style  <<<<<<<\n");
		  ts << tr("*>>>>>>>   for ") << molBaseName;
		  ts << tr(" which generated by Avogadro.          <<<<<<<\n\n");
		  ts << "BONDS\n";
		  ts << "!\n";
          ts << "!V(bond) = Kb(b - b0)**2\n";
          ts << "!\n";
          ts << "!Kb: kcal/mole/A**2\n";
          ts << "!b0: A\n";
          ts << "!\n";
          ts << "!atom type Kb          b0\n";
          ts << "!\n";
		 
		 for( int i = 0; i < bkmol.size(); i++ )
	     {

		  ts << qSetRealNumberPrecision(5) << fixed << right
			 << bkmol[i].aT1 << " "
             <<  bkmol[i].aT2 << "    "
			 << qSetRealNumberPrecision(9)
			 <<  bkmol[i].Kb << " "
			 <<  bkmol[i].req << endl;
   	     }
		 
		  ts << endl; // insert empty line

		 ts << "ANGLES\n";
		 ts << "!\n";
		 ts << "!V(angle) = Ktheta(Theta - Theta0)**2\n";
		 ts << "!\n";
		 ts << "!V(Urey-Bradley) = Kub(S - S0)**2\n";
		 ts << "!\n";
		 ts << "!Ktheta: kcal/mole/rad**2\n";
		 ts << "!Theta0: degrees\n";
		 ts << "!Kub: kcal/mole/A**2 (Urey-Bradley)\n";
		 ts << "!S0: A\n";
		 ts << "!\n";
		 ts << tr("!atom types     Ktheta    Theta0   Kub     S0\n");
		 ts << "!\n";


		 for( int i = 0; i < akmol.size(); i++ )
	     {
		    ts  << qSetRealNumberPrecision(5) << fixed << right
				<<  akmol[i].aT1 << " "
               <<  akmol[i].aT2 << " "
			   <<  akmol[i].aT3 << "    "
			   <<  qSetRealNumberPrecision(9)
			   <<  akmol[i].Kth << " "
			   <<  akmol[i].theq << endl;
	     
		 }


		ts << endl; // an empty line

	    ts << "DIHEDRALS" << endl;
        ts << "!\n";
        ts << "!V(dihedral) = Kchi(1 + cos(n(chi) - delta))\n";
        ts << "!\n";
        ts << "!Kchi: kcal/mole\n";
        ts << "!n: multiplicity\n";
        ts << "!delta: degrees\n";
        ts << "!\n";
        ts << "!atom types             Kchi    n   delta\n";
		ts << "!\n";


	   for( int i = 0; i < tkmol.size(); i++ )
	   {
		      ts  << qSetRealNumberPrecision(5) << fixed << right
			      <<  tkmol[i].aT1 << " "
                  <<  tkmol[i].aT2 << " "
			      <<  tkmol[i].aT3 << " "
                  <<  tkmol[i].aT4 << "    "
			      <<  tkmol[i].vn2 << " "
			      <<  qSetRealNumberPrecision(0) 
			      <<  tkmol[i].n << "  "
			      <<  qSetRealNumberPrecision(5) << fixed << right
			      <<  tkmol[i].gamma << endl;
   	    }


	    ts << endl; // an empty line

	    ts << "IMPROPER" << endl;
        ts << "!V(improper) = Kpsi(psi - psi0)**2\n";
        ts << "!\n";
        ts << "!Kpsi: kcal/mole/rad**2\n";
        ts << "!psi0: degrees\n";
        ts << "!note that the second column of numbers (0) is ignored\n";
        ts << "!\n";

      /*
	   for( int i = 0; i < okmol.size(); i++ )
	  {
		    ts << qSetRealNumberPrecision(5) << fixed << right
		    	<<  okmol[i].aT1 << " "
               <<  okmol[i].aT2 << " "
			   <<  okmol[i].aT3 << " "
               <<  okmol[i].aT4 << "    "
			   <<  okmol[i].vn2 << "   "
			   << "0    "
			   <<  okmol[i].gamma << endl;

   	  }
      */

	// ts << endl; // an empty line

	ts << "NONBONDED" << endl;
	ts << "!\n";
	ts << "!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]\n";
	ts << "!\n";
	ts << "!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)\n";
	ts << "!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j\n";
	ts << "!\n";
	ts << "!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4\n";
	ts << "!\n";


	for( int i = 0; i< nkmol.size() ; i++ )
	{

		ts  << qSetRealNumberPrecision(5) << fixed << right
			<< nkmol[i].aT << "     "
			<< "0.00000" << "   "
			<< nkmol[i].eps <<  "   "
			<< nkmol[i].R
			<< endl;
	}



	ts << endl; // an empty line

	ts << "END" << endl << endl; 

	  }// end( if(m_molecule != NULL)

	  else
	  {
		  ts << tr("! No molecule was found.\n");
		  ts << tr("! These psf and pot preview panes are NOT automatically updated !\n");
		  ts << tr("! Please click the update button to update !\n");
	  }

	  return buffer;

   }// end potGhemical()


	void PsfPotDialog::setForceField(int n)
	{
		switch(n)
		{
		  case 0:
          default:
			  m_forceField = "GAFF";
			  break;

		  case 1:
			  m_forceField = "UFF";
			  break;

		  case 2:
			  m_forceField = "Ghemical";
			  break;
		}

		updatePreviewText();

	}


}