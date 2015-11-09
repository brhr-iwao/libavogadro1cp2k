/**********************************************************************
  PsfPotDialog - Charmm-style Psf/Pot Generator Dialog

  Copyright (C) aoyama iwao

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

#include "psfpotdialog.h"

#include <avogadro/molecule.h>
#include <avogadro/atom.h>

//#include <openbabel/mol.h>

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
//#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>


#include <QString>
#include <QTextStream>
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>

using namespace OpenBabel;
using namespace std;

namespace Avogadro
{
  PsfPotDialog::PsfPotDialog(QWidget *parent, Qt::WindowFlags f) : QDialog(parent, f), m_molecule(0), m_savePath("")
  {
	  ui.setupUi(this);

	  connect(ui.closeButton, SIGNAL(clicked()), this, SLOT(close()) );
	  connect(ui.psfGenerateButton, SIGNAL(clicked()), this, SLOT( psfGenerateClicked() ) );
	  connect(ui.potGenerateButton, SIGNAL(clicked()), this, SLOT( potGenerateClicked() ) );
	  connect(ui.updatePreviewButton, SIGNAL(clicked()), this, SLOT( updatePreviewText() ) );

	  Updated = true;

	  updatePreviewText();
  }

  PsfPotDialog::~PsfPotDialog()
  {
  }

  void PsfPotDialog::setMolecule(Molecule *molecule)
  {
	  // Disconnect the old molecule first...
	  if (m_molecule) disconnect(m_molecule, 0, this, 0);

	  Updated = false;
   
       m_molecule = molecule;
  }

  void PsfPotDialog::readSettings(QSettings&)
  {
  }

  void PsfPotDialog::writeSettings(QSettings&) const
  {
  }

  QString PsfPotDialog::psfPreviewPane()
  {
	  return psfGaff();
  }

  QString PsfPotDialog::potPreviewPane()
  {
	  return potGaff();
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
          msgBox.setInformativeText(tr("Do you want to update before generating file?"));
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
                          tr("Save Psf File"), QString("psf"));
  
  }

  void PsfPotDialog::potGenerateClicked()
  {
	  if(!Updated)
	  {
		  QMessageBox msgBox;
          msgBox.setText(tr("The moleule has been modified."));
          msgBox.setInformativeText(tr("Do you want to update before generating file?"));
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
                          tr("Save Pot File"), QString("pot"));

  
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

  // convert a first and a second letter of a character array
  // into the integer-type acsii code
  int PsfPotDialog::chrs2int( const char* chrs )
  {
	int retval = 0;

    char fstchr[32], scdchr[32];
	int fstval, scdval; 

	sprintf( fstchr, "%i", chrs[0] );
    sprintf( scdchr, "%i", chrs[1] ); // chrs[1]=='\0' case may occur.

	fstval = atoi(fstchr);
	scdval = atoi(scdchr);

	retval = fstval*100 + scdval;

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


  bool PsfPotDialog::torTypeEqual( const torKind& left, const torKind& right )
  {
	if( ( strcmp( left.aT1, right.aT1 ) == 0 ) && (strcmp( left.aT2, right.aT2 ) == 0)
		&& (strcmp( left.aT3, right.aT3 ) == 0) && (strcmp( left.aT4, right.aT4 ) == 0) )
		return true;

	// reverse order (1-2-3-4 to 4-3-2-1)
	else if( ( strcmp( left.aT1, right.aT4 ) == 0 ) && (strcmp( left.aT2, right.aT3 ) == 0)
		&& (strcmp( left.aT3, right.aT2 ) == 0) && (strcmp( left.aT4, right.aT1 ) == 0) )
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


  /*
  // Calulate bond angle force constant for averaged distances and angles
  // of non-listed triple atom kinds in input molecule.
  // See J.Wang et al, J.Comp.Chem., 25, 1157-1174(2004)
  */
  void PsfPotDialog::acquireWholeAngles( OBMol mol, vector<angleKind>* ak )
  {
	ak->clear();
	
   FOR_ANGLES_OF_MOL(angle, mol)
	{
		OBAtom* a = mol.GetAtom((*angle)[1] + 1);
		OBAtom* b = mol.GetAtom((*angle)[0] + 1); // vertex
		OBAtom* c = mol.GetAtom((*angle)[2] + 1);

	    OBPairData* aaType = (OBPairData*) a->GetData("FFAtomType");
		OBPairData* baType = (OBPairData*) b->GetData("FFAtomType");
		OBPairData* caType = (OBPairData*) c->GetData("FFAtomType");

		angleKind tempAngle;

		if( chrs2int(aaType->GetValue().c_str()) <= 
              chrs2int(caType->GetValue().c_str() ) )
		{
		   strcpy( tempAngle.aT1, aaType->GetValue().c_str());
		   strcpy( tempAngle.aT2, baType->GetValue().c_str());
		   strcpy( tempAngle.aT3, caType->GetValue().c_str());
		}

		else
		{
		   strcpy( tempAngle.aT1, caType->GetValue().c_str());
		   strcpy( tempAngle.aT2, baType->GetValue().c_str());
		   strcpy( tempAngle.aT3, aaType->GetValue().c_str());
		}

		tempAngle.Kth = 0.0;
		tempAngle.theq = 0.0;

		if( a->GetAtomicNum() <= c-> GetAtomicNum() ) // store ascending order
		{
			tempAngle.aN1 = a->GetAtomicNum();
			tempAngle.aN2 = b->GetAtomicNum();
			tempAngle.aN3 = c->GetAtomicNum();
		}

		else
		{
			tempAngle.aN1 = c->GetAtomicNum();
			tempAngle.aN2 = b->GetAtomicNum();
			tempAngle.aN3 = a->GetAtomicNum();
		}

		tempAngle.r1c = a->GetDistance( b );
		tempAngle.r2c = b->GetDistance( c );
        tempAngle.thc = b->GetAngle( a, c );

		ak->push_back( tempAngle );
	}

	return;
  }

  /*
  // Calulate bond angle force constant for averaged distances and angles
  // of non-listed triple atom kinds in input molecule.
  // See J.Wang et al, J.Comp.Chem., 25, 1157-1174(2004)
  */
  double PsfPotDialog::gaffAngleC( unsigned int atomicNum )
  {
	switch( atomicNum )
	{
	  case 6: // carbon
		return 1.339;
	  case 7: // nitrogen
		  return 1.300;
	  case 8: // oxygen
		  return 1.249;
	  case 15: // phosphorus
		  return 0.906;
	  case 16: // sulfur
		  return 1.448;
	  default:
		  return 0.0;
	}
  }

   /*
  // Calulate bond angle force constant for averaged distances and angles
  // of non-listed triple atom kinds in input molecule.
  // See J.Wang et al, J.Comp.Chem., 25, 1157-1174(2004)
  */
  double PsfPotDialog::gaffAngleZ( unsigned int atomicNum )
  {
	switch( atomicNum )
	{
	  case 1: // hydrogen
		  return 0.784;
	  case 6: // carbon
		  return 1.183;
	  case 7: // nitrogen
		  return 1.212;
	  case 8: // oxygen
		  return 1.219;
	  case 9: // fluorine
		  return 1.166;
	  case 17: // chlorine
		  return 1.272;
	  case 35: // bromine
		  return 1.378;
	  case 53: // iodine
		  return 1.398;
	  case 15: // phosphorus
		  return 1.620;
	  case 16: // sulfur
		  return 1.280;
	  default:
		  return 0.0;
	}
  }

  /*
  // The following functions will be move to the other cpp files
  //
  //
  //
  */

  	QString PsfPotDialog::psfGaff()
	{
	  QString buffer;
      QTextStream ts(&buffer);

	  QString molBaseName;
	  OBMol mol;

	  QString debugInfo;
	  QTextStream dits(&debugInfo);// debugInfo text stream;

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

	       OBForceField* pFF = OBForceField::FindForceField("GAFF");
           if (!pFF)
	       {
               // qDebug() << "Could not find GAFF!"; 
               return tr("! Could not find GAFF!");
	       }

            if (!pFF->Setup(mol))
	       {
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
		    // mols[molid].AddHydrogens(false, true); // hydrogens must be added before Setup(mol) is called
		    // after calling AddHydrogen(), atom indices are permutated.

	           if (!pFF->Setup( mols[molid] ))
	           {
                     dits << tr("! Could not setup GAFF for the fragment of ") << molid << ".\n";
		             dits << tr("! It may go well if you add hydrogens to the fragment of ") << molid << ".\n";
                     return debugInfo;
			   }

			   if(!pFF->GetAtomTypes(mols[molid] ))
              {
	               dits << tr(" could not find atom types of the fragment of") << molid << ".\n";
	               return debugInfo;
               }
            }

			ts << "PSF\n\n";
            ts << "      3 !NTITLE\n";
			ts << "      " << left << molBaseName << ".\n";
            ts << tr("      Avogadro generated Protein Structure File (PSF)\n");
			ts << tr("      by using General AMBER Force Fields (GAFF).\n\n");

			// ts << qSetFieldWidth(7) << left << mol.NumAtoms() << "\n"<< qSetFieldWidth(0);
			ts << qSetFieldWidth(5) << right << mol.NumAtoms() << " !NATOMS\n";

   	      int acum_natom = 0;
	      int acum_nres = 0;

           for ( molid = 0; molid < mols.size() ; molid++ )
	      {
		      FOR_ATOMS_OF_MOL( atom, mols[molid] )
		     {
			      int atomid = atom->GetId() + acum_natom; // GetId() retrives unique atom index of mol and mols starts from 0. 
			                                               // (GetId() was different from GetIdx() 
			                                               // which starts from 1 and are changable by AddHyrogen() or Separate())

			      OBPairData *type = (OBPairData*) atom->GetData("FFAtomType");

				  ts << qSetFieldWidth(7) << right << atomid + 1 << qSetFieldWidth(0)
					 << " SEG" << left << molid + 1;
				  ts << qSetFieldWidth(7) << right 
					  << atom->GetResidue()->GetIdx() + 1 + acum_nres;
				  ts << qSetFieldWidth(5) << right << atom->GetResidue()->GetName().c_str() 
					  << " " << qSetFieldWidth(0);
				  ts << left << etab.GetSymbol(atom->GetAtomicNum()) << "  ";
				  ts << left << type->GetValue().c_str() << "    ";
				  ts << qSetRealNumberPrecision(5) << right << atom->GetPartialCharge() <<  "  ";
				  ts << qSetRealNumberPrecision(5) << right << atom->GetExactMass() << "      0\n";

			 /*
			 printf("%7d SEG%d", atomid + 1, molid + 1 );
			 printf("%7d", atom->GetResidue()->GetIdx() + 1 + acum_nres );
			 printf("%5s ", atom->GetResidue()->GetName().c_str() );
	         printf("%-s ", etab.GetSymbol(atom->GetAtomicNum()));
			 printf("%-s    ", type->GetValue().c_str() ); // x4
	         printf("%.5f  ", atom->GetPartialCharge());
	         printf("%.5f      0\n", atom->GetExactMass() );
			 */

		      }	   
		          acum_natom += mols[molid].NumAtoms();
		          acum_nres += mols[molid].NumResidues();
	      }

	      ts << endl;
          ts << qSetFieldWidth(5) << right << mol.NumBonds() << " !NBONDS\n";
		  // printf("%5d !NBOND bonds\n", mol.NumBonds() );

		 // bond pairs
         int i = 0;
         acum_natom = 0;
         for ( molid = 0; molid < mols.size() ; molid++ )
         {
	          FOR_BONDS_OF_MOL(bond, mols[molid] )
             {
	              ts << "      " 
				     << bond->GetBeginAtom()->GetId() + 1 + acum_natom 
				     << "      " 
				    << bond->GetEndAtom()->GetId() + 1 + acum_natom ;
	               i++;

	                if( i >= 4 ) // four pairs of atoms per line
	              {
			             ts << endl;
	                     i=0;
	              }
	              else continue;
              }

		         acum_natom += mols[molid].NumAtoms();
          }
          // end of bonds

		  ts << endl << endl;

		  // number of angles
          int nAngle = 0;
          FOR_ANGLES_OF_MOL(angle, mol){ nAngle++; }

          ts << qSetFieldWidth(5) << right << nAngle << " !NTHATA: angles" << endl;

		  // angle triples
          OBAtom *a, *b, *c;
           i = 0;
           acum_natom = 0;
           for ( molid = 0; molid < mols.size() ; molid++ )
           {
	            FOR_ANGLES_OF_MOL( angle, mols[molid] )
		        {
		             a = mols[molid].GetAtom((*angle)[1] + 1);
	                 b = mols[molid].GetAtom((*angle)[0] + 1); // vertex
	                 c = mols[molid].GetAtom((*angle)[2] + 1);

	                ts << "      " 
			        	<< a->GetId() + 1 + acum_natom
			    	    << "      " 
				        << b->GetId() + 1 + acum_natom
				         << "      " 
				       << c->GetId() + 1 + acum_natom;
                     i++;

	                if( i >= 3 ) // three triples of atoms per line
	                { 
		                ts << endl;
		                 i = 0;
	                }
	                else continue;

		           }

	           acum_natom += mols[molid].NumAtoms();
              }

          ts << endl;
         // end of angles

		  ts << endl; // an empty line

        // number of torsion (dihedral) quadreples
        int nTorsion = 0;
        FOR_TORSIONS_OF_MOL(torsion, mol){ nTorsion++ ;}
		ts << qSetFieldWidth(5) << right << nTorsion << " !NPHI: dihedrals\n";

        // torsion (dihedral) quadreples
        OBAtom* d; // OBAtom *a, *b and *c are already decleard above.
        i = 0;
        acum_natom = 0;
        for ( molid = 0; molid < mols.size() ; molid++ )
        {
	         FOR_TORSIONS_OF_MOL(torsion, mols[molid] )
           {
	           a = mols[molid].GetAtom((*torsion)[0] + 1);
	           b = mols[molid].GetAtom((*torsion)[1] + 1);
	           c = mols[molid].GetAtom((*torsion)[2] + 1);
	           d = mols[molid].GetAtom((*torsion)[3] + 1);

	           ts << "      " 
				  << a->GetId() + 1 + acum_natom
				  << "      " 
				  << b->GetId() + 1 + acum_natom
				  << "      "
				  << c->GetId() + 1 + acum_natom 
				  << "      " <<
				     d->GetId() + 1 + acum_natom ;
	           i++;

	          if( i >= 2 ) // two quaruples of atoms per line
	          { 
		        ts << endl;
		        i = 0;
	           }
	        else continue;
           }

           acum_natom += mols[molid].NumAtoms();
         }
         // end of dihedrals

          ts << endl; // an empty line
 
          // number of out-of-plane ( improper ) dihedral quadreples
          int nOOP = 0;
          FOR_ATOMS_OF_MOL (atom, mol) 
         { 
	        if( atom->GetImplicitValence() == 3 ) {nOOP++;}
         }
         ts << nOOP << " !NIMPHI: impropers" << endl;
         // printf("%5d !NIMPHI: impropers\n", nOOP);

        // out-of-plane (improper) dihedral quareples
         i = 0;
        acum_natom = 0;
       for ( molid = 0; molid < mols.size() ; molid++ )
       {
	       FOR_ATOMS_OF_MOL ( atom, mols[molid] ) 
           {
	          if( atom->GetImplicitValence() == 3 )
	          {
		            ts << "      " << atom->GetId() + 1 + acum_natom;

                   FOR_NBORS_OF_ATOM( nbr, &*atom )
		          {

			         ts << "      " << nbr->GetId() + 1 + acum_natom;
		           }
		          i++;

		          if( i >= 2 ) // two quaruples of atoms per line
		    	  {
		             ts << endl;
				     i = 0;
		          }

			  }
	           else continue;
		    }

		acum_natom += mols[molid].NumAtoms();
       }
       // end of out-of-plane

        ts << endl; // an empty line

        ts << "    0 !NDON\n\n";
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

	  QString debugInfo;
	  QTextStream dits(&debugInfo); // debugInfo text stream;

	  int i;

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

		   // Acquire bond kind
           vector<bondKind> bkmol;

           FOR_BONDS_OF_MOL(bond, mol)
	       {
	       	 OBAtom* beginAtom = bond->GetBeginAtom();
	    	 OBAtom* endAtom = bond->GetEndAtom();

	         OBPairData *beginType = (OBPairData*) beginAtom->GetData("FFAtomType");
		     OBPairData *endType = (OBPairData*) endAtom->GetData("FFAtomType");

	         bondKind tempKind;

		       // order two atom type members as ascening.
		       if( chrs2int(beginType->GetValue().c_str()) <= 
                 chrs2int(endType->GetValue().c_str() ) )
		        {
		           strcpy( tempKind.aT1, beginType->GetValue().c_str());
		           strcpy( tempKind.aT2, endType->GetValue().c_str());
		        }

		       else
		       {
		          strcpy( tempKind.aT1, endType->GetValue().c_str());
		          strcpy( tempKind.aT2, beginType->GetValue().c_str());
		       }

		      tempKind.Kb = 0.0;
		      tempKind.req = 0.0;

		      bkmol.push_back( tempKind );
		
	       }

		   sort( bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeComp );
	       vector<bondKind>::iterator bondEnd = unique(bkmol.begin(), bkmol.end(), Avogadro::PsfPotDialog::bondTypeEqual );
	       bkmol.erase( bondEnd, bkmol.end() );

		   	// acquire angle kind
	        vector<angleKind> akmol;

            acquireWholeAngles( mol, &akmol );
            sort( akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeComp );
	        vector<angleKind>::iterator angleEnd = unique(akmol.begin(), akmol.end(), Avogadro::PsfPotDialog::angleTypeEqual );
	        akmol.erase( angleEnd, akmol.end() );
			// end (acquire angle kind)

	        //  aquire torsion kind
	        vector<torKind> tkmol;

            FOR_TORSIONS_OF_MOL(torsion, mol)
            {
           	    OBAtom* a = mol.GetAtom((*torsion)[0] + 1);
	            OBAtom* b = mol.GetAtom((*torsion)[1] + 1);
	            OBAtom* c = mol.GetAtom((*torsion)[2] + 1);
	            OBAtom* d = mol.GetAtom((*torsion)[3] + 1);

	            OBPairData* atType = (OBPairData*) a->GetData("FFAtomType");
	            OBPairData* btType = (OBPairData*) b->GetData("FFAtomType");
	            OBPairData* ctType = (OBPairData*) c->GetData("FFAtomType");
	            OBPairData* dtType = (OBPairData*) d->GetData("FFAtomType");

	             torKind tempTor;

	            if( chrs2int(atType->GetValue().c_str()) <= 
                      chrs2int(dtType->GetValue().c_str() ) )
	             {
		             strcpy( tempTor.aT1, atType->GetValue().c_str());
		             strcpy( tempTor.aT2, btType->GetValue().c_str());
		             strcpy( tempTor.aT3, ctType->GetValue().c_str());
		             strcpy( tempTor.aT4, dtType->GetValue().c_str());
	             }

	            else
	            {
		             strcpy( tempTor.aT1, dtType->GetValue().c_str());
		             strcpy( tempTor.aT2, ctType->GetValue().c_str());
		             strcpy( tempTor.aT3, btType->GetValue().c_str());
		             strcpy( tempTor.aT4, atType->GetValue().c_str());
	             }

	            tempTor.m = 0;
	            tempTor.vn2 = 0.0;
	            tempTor.gamma = 0.0;
	            tempTor.n = 0.0;

	            tkmol.push_back( tempTor );

              }

	          sort( tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeComp );
	          vector<torKind>::iterator torEnd = unique(tkmol.begin(), tkmol.end(), Avogadro::PsfPotDialog::torTypeEqual );
	          tkmol.erase( torEnd, tkmol.end() );
	          // end (aqire torsion kind)

			  // acquire oopKind (out-of-plane improper dihedral )
	         vector<oopKind> okmol;

             FOR_ATOMS_OF_MOL (atom, mol) 
             {
	            if( atom->GetImplicitValence() == 3 )
	            {
		            OBAtom* pa[4];

		            pa[0] = (OBAtom*)&*atom;

                    i=1;
		            FOR_NBORS_OF_ATOM( nbr, &*atom )
		              {  pa[i] = (OBAtom*)&*nbr; i++; }

	                OBPairData* opa1Type = (OBPairData*) pa[0]->GetData("FFAtomType");
	                OBPairData* opa2Type = (OBPairData*) pa[1]->GetData("FFAtomType");
	                OBPairData* opa3Type = (OBPairData*) pa[2]->GetData("FFAtomType");
	                OBPairData* opa4Type = (OBPairData*) pa[3]->GetData("FFAtomType");

		            oopKind tempOop;

		            strcpy( tempOop.aT1, opa1Type->GetValue().c_str());

		            if(  ( chrs2int(opa2Type->GetValue().c_str() ) <= chrs2int(opa4Type->GetValue().c_str() ) ) &&
                         ( chrs2int(opa4Type->GetValue().c_str() ) <= chrs2int(opa3Type->GetValue().c_str() ) ) )
		            {
		                 strcpy( tempOop.aT2, opa2Type->GetValue().c_str());
				         strcpy( tempOop.aT3, opa4Type->GetValue().c_str());
				         strcpy( tempOop.aT4, opa3Type->GetValue().c_str());
		             }

		             else if( ( chrs2int(opa3Type->GetValue().c_str()) <= chrs2int(opa2Type->GetValue().c_str() ) ) &&
			                  ( chrs2int(opa2Type->GetValue().c_str() ) <= chrs2int(opa4Type->GetValue().c_str() ) ) )
		             {
		                  strcpy( tempOop.aT2, opa3Type->GetValue().c_str());
				          strcpy( tempOop.aT3, opa2Type->GetValue().c_str());
				          strcpy( tempOop.aT4, opa4Type->GetValue().c_str());
		             }

		             else if( ( chrs2int(opa4Type->GetValue().c_str()) <= chrs2int(opa2Type->GetValue().c_str() ) ) &&
                          ( chrs2int(opa2Type->GetValue().c_str() ) <= chrs2int(opa3Type->GetValue().c_str() ) ) )
		             {
		                     strcpy( tempOop.aT2, opa4Type->GetValue().c_str());
				             strcpy( tempOop.aT3, opa2Type->GetValue().c_str());
				             strcpy( tempOop.aT4, opa3Type->GetValue().c_str());
		              }

		             else if( ( chrs2int(opa3Type->GetValue().c_str()) <= chrs2int(opa2Type->GetValue().c_str() ) ) &&
                           ( chrs2int(opa2Type->GetValue().c_str() ) <= chrs2int(opa4Type->GetValue().c_str() ) ) )
		              {
		                     strcpy( tempOop.aT2, opa3Type->GetValue().c_str());
				             strcpy( tempOop.aT3, opa2Type->GetValue().c_str());
				             strcpy( tempOop.aT4, opa4Type->GetValue().c_str());
		              }

		            else if( ( chrs2int(opa4Type->GetValue().c_str()) <= chrs2int(opa3Type->GetValue().c_str() ) ) &&
                          ( chrs2int(opa3Type->GetValue().c_str() ) <= chrs2int(opa2Type->GetValue().c_str() ) ) )
		             {
		                     strcpy( tempOop.aT2, opa4Type->GetValue().c_str());
				             strcpy( tempOop.aT3, opa3Type->GetValue().c_str());
				             strcpy( tempOop.aT4, opa2Type->GetValue().c_str());
		              }

		             else
		             {
			                 strcpy( tempOop.aT2, opa2Type->GetValue().c_str());
				             strcpy( tempOop.aT3, opa3Type->GetValue().c_str());
				             strcpy( tempOop.aT4, opa4Type->GetValue().c_str());
		              }

		            tempOop.vn2 = 0.0;
		            tempOop.gamma = 0.0;
		            tempOop.n = 0;

		           okmol.push_back( tempOop );

	             }

	           else continue;
             }

	         sort( okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeComp );
	         vector<oopKind>::iterator oopEnd = unique(okmol.begin(), okmol.end(), Avogadro::PsfPotDialog::oopTypeEqual );
	         okmol.erase( oopEnd, okmol.end() );
              // end (acquire oopkind ) (out-of-plane improper dihedral )

			  // acquire nbkind (non boning interaction)
              vector<nbKind> nkmol;

	          FOR_ATOMS_OF_MOL (atom, mol)
	          {
	                  OBPairData *tp = (OBPairData*) atom->GetData("FFAtomType");

		              nbKind tempNb;

		              strcpy( tempNb.aT, tp->GetValue().c_str() );
		              tempNb.R = 0.0;
		              tempNb.eps = 0.0;

		              nkmol.push_back( tempNb );
	           }

	          sort( nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeComp );
	          vector<nbKind>::iterator nbEnd = unique(nkmol.begin(), nkmol.end(), Avogadro::PsfPotDialog::nbTypeEqual );
	          nkmol.erase( nbEnd, nkmol.end() );
	          // end (acquire nbkind) (non boning interaction)

			  // Acquire gaff.dat file data
	          char f[] = "gaff.dat";
	          fstream s(f, ios::in);

	          if(s.fail())
				  return "Can't open gaff.dat";

	          char buffer[BUFF_SIZE];
	          char* pchar;

			  int iBlock = 0; // 0: atom block, 1: bond block, 
	                            // 2: angle block 3: torsion block, 4: improper block
	          bool IsNbBlock = false; // Is in non-bonding block currently?

	          while (1) // read a file...
             {
                  if (!s.getline(buffer,BUFF_SIZE)) // file f ends.
                       break;

		          // if an empty line (return only line) is detected, increment iBlock.
		          if( strcmp(buffer, "") == 0 ) 
			          iBlock++; 

		          if( strcmp(buffer, "MOD4      RE") == 0 ) 
			            IsNbBlock = true; // if start of non-bonding block is detected, switch the flag true

		          if( strcmp(buffer, "END") == 0 )
			             IsNbBlock = false; // if end of non-bonding block is detected, switch the flag false

		         // read bond block
                  if( iBlock == 1 )
		          {
			           char tAtom1[3], tAtom2[3]; // temporarily stored atom types read from lines.
			           strcpy( tAtom1, "" ); // initialize
			           strcpy( tAtom2, "" ); // initialize

			           double tKb = 0.0;
			           double tReq = 0.0;

		               pchar = strtok( buffer, "-, "); // call first token ( separated by a hyphen or blank ) (first atom type)
			           if( pchar != NULL )
			               { strcpy( tAtom1, pchar );}
  
			           pchar = strtok( NULL, "-, "); // call second token (second atom typem)
			           if( pchar != NULL )
			               { strcpy( tAtom2, pchar );}

			           pchar = strtok( NULL, "-, "); // call third token (force constant)
			          if( pchar != NULL)
			                 { tKb = atof(pchar);}

			           pchar = strtok( NULL, "-, "); // call fourth token (equibrium bond length)
			           if( pchar != NULL )
			           { tReq = atof(pchar);}

			          // if match temporary atom types and mol bonds ones
			          for( i = 0 ; i < bkmol.size() ; i++ )
			          {
				            if ( ( stricmp( bkmol[i].aT1, tAtom1 ) ==0 && stricmp( bkmol[i].aT2, tAtom2) == 0 ) ||
                                 ( stricmp( bkmol[i].aT1, tAtom2 ) ==0 && stricmp( bkmol[i].aT2, tAtom1) == 0 ) )
				            {
					              bkmol[i].Kb = tKb;
					              bkmol[i].req = tReq;
				             }
			           }

		          }
		         // end ( read bond block )

	            // read angle block
               else if( iBlock == 2 )
	          {
		   	       char tAtom1[3], tAtom2[3], tAtom3[3]; // temporarily stored atom types read from lines.
			       strcpy( tAtom1, "" ); // initialize
			       strcpy( tAtom2, "" ); // initialize
			       strcpy( tAtom3, "" ); // initialize

			       double tKth = 0.0;
			       double tTheq = 0.0;

			       pchar = strtok( buffer, "-, "); // call first token ( separated by a hyphen or blank ) (first atom type)
			       if( pchar != NULL )
			          { strcpy( tAtom1, pchar );}
  
			       pchar = strtok( NULL, "-, "); // call second token (second atom type)
			       if( pchar != NULL )
			            { strcpy( tAtom2, pchar );}

			       pchar = strtok( NULL, "-, "); // call third token (third atom type)
			       if( pchar != NULL )
			           { strcpy( tAtom3, pchar );}

			       pchar = strtok( NULL, "-, "); // call fourth token (bending constant)
			       if( pchar != NULL)
			           { tKth = atof(pchar);}

			       pchar = strtok( NULL, "-, "); // call fifth token (equibrium bend angle)
			       if( pchar != NULL )
			           { tTheq = atof(pchar);}

			       for( i = 0 ; i < akmol.size() ; i++ )
			       {
				     if ( ( stricmp( akmol[i].aT1, tAtom1 ) ==0 && stricmp( akmol[i].aT2, tAtom2) ==0  && stricmp( akmol[i].aT3, tAtom3) == 0 ) ||
                       ( stricmp( akmol[i].aT1, tAtom3 ) ==0 && stricmp( akmol[i].aT2, tAtom2) ==0  && stricmp( akmol[i].aT3, tAtom1) == 0 )  )
				     {
					     akmol[i].Kth = tKth;
					     akmol[i].theq = tTheq;
				      }
			        }
	          }
	         // end ( read angle block)

	         // read torsion dihedral block
	         else if( iBlock == 3 )
	        {
		        char tAtom1[3], tAtom2[3], tAtom3[3], tAtom4[3]; // temporarily stored atom types read from lines.
			    strcpy( tAtom1, "" ); // initialize
			    strcpy( tAtom2, "" ); // initialize
			    strcpy( tAtom3, "" ); // initialize
			    strcpy( tAtom4, "" ); // initialize

	            double  tVn2 = 0.0; 
	            double  tGamma = 0.0; 
                double  tN = 0.0;

			    pchar = strtok( buffer, "-, "); // call first token ( separated by a hyphen or blank ) (first atom type)
			    if( pchar != NULL )
			       { strcpy( tAtom1, pchar );}

			    pchar = strtok( NULL, "-, "); // call second token ( second atom type )
			    if( pchar != NULL )
			       { strcpy( tAtom2, pchar );}

			    pchar = strtok( NULL, "-, "); // call third token( third atom type )
			    if( pchar != NULL )
			       { strcpy( tAtom3, pchar );}

			    pchar = strtok( NULL, "-, "); // call fourth token ( fourth atom type )
			    if( pchar != NULL )
			       { strcpy( tAtom4, pchar );}

			    pchar = strtok( NULL, "-, "); // ignore fifth column (number of path)...
			    pchar = strtok( NULL, "-, "); // and call sixth token ( magnitude of torsion )
			    if( pchar != NULL)
			       { tVn2 = atof(pchar);} // magnitude of torsion in kcal/mol (Vn/2)

			    pchar = strtok( NULL, "-, "); // call seventh token ( phase offset )
			    if( pchar != NULL)
			        { tGamma = atof(pchar);} // phese offset in degree 

			    pchar = strtok( NULL, "-, "); // call eightth token ( periodicity of torsion )
			    if( pchar != NULL)
			       { tN = atof(pchar);} // periodicity of torsion

			    for( i = 0 ; i < tkmol.size() ; i++ )
			    {
                   // temporarily only 2-3 atoms are considered
				  if ( ( stricmp( tkmol[i].aT2, tAtom2 ) ==0 && stricmp( tkmol[i].aT3, tAtom3) == 0 ) ||
                        ( stricmp( tkmol[i].aT2, tAtom3 ) ==0 && stricmp( tkmol[i].aT3, tAtom2) == 0 ) )
				   {
					 tkmol[i].vn2 = tVn2;
					 tkmol[i].gamma = tGamma;
					 tkmol[i].n = tN;
				   }
			    }

	          }
           // end ( read torsion dihedral block ) 

          // read improper dihedral block
	      else if( iBlock == 4 )
	      {
		      char tAtom1[3], tAtom2[3], tAtom3[3], tAtom4[3]; // temporarily stored atom types read from lines.
			  strcpy( tAtom1, "" ); // initialize
			  strcpy( tAtom2, "" ); // initialize
			  strcpy( tAtom3, "" ); // initialize
			  strcpy( tAtom4, "" ); // initialize

			  double  tVn2 = 0.0; 
	          double  tGamma = 0.0; 
              double  tN = 0.0;

			  pchar = strtok( buffer, "-, "); // call first token ( separated by a hyphen or blank ) (first atom type)
			  if( pchar != NULL )
			  { strcpy( tAtom1, pchar );}

			  pchar = strtok( NULL, "-, "); // call second token ( second atom type )
			  if( pchar != NULL )
			  { strcpy( tAtom2, pchar );}

			   pchar = strtok( NULL, "-, "); // call third token( third atom type )
			  if( pchar != NULL )
			  { strcpy( tAtom3, pchar );}

			  pchar = strtok( NULL, "-, "); // call fourth token ( fourth atom type )
			  if( pchar != NULL )
			  { strcpy( tAtom4, pchar );}

			  pchar = strtok( NULL, "-, "); // and call fifth token ( magnitude of torsion )
			  if( pchar != NULL)
			  { tVn2 = atof(pchar);} // magnitude of torsion in kcal/mol (Vn/2)

			  pchar = strtok( NULL, "-, "); // call sixth token ( phase offset )
			  if( pchar != NULL)
			  { tGamma = atof(pchar);} // phese offset in degree 

			  pchar = strtok( NULL, "-, "); // call seventh token ( periodicity of torsion )
			  if( pchar != NULL)
			  { tN = atof(pchar);} // periodicity of torsion

			   for( i = 0 ; i < okmol.size() ; i++ )
			   {
				  if( stricmp( okmol[i].aT1, tAtom1 ) ==0 || stricmp(tAtom1, "X" ) == 0 )
				  {
					  if( ( stricmp( okmol[i].aT2, tAtom2 ) == 0 && stricmp( okmol[i].aT3, tAtom3 ) == 0 && stricmp( okmol[i].aT4, tAtom4 ) == 0 )
						  || ( stricmp( okmol[i].aT2, tAtom2 ) == 0 && stricmp( okmol[i].aT3, tAtom4 ) == 0 && stricmp( okmol[i].aT4, tAtom3 ) == 0 )
						  || ( stricmp( okmol[i].aT2, tAtom3 ) == 0 && stricmp( okmol[i].aT3, tAtom2 ) == 0 && stricmp( okmol[i].aT4, tAtom4 ) == 0 )
						  || ( stricmp( okmol[i].aT2, tAtom3 ) == 0 && stricmp( okmol[i].aT3, tAtom4 ) == 0 && stricmp( okmol[i].aT4, tAtom2 ) == 0 )
						  || ( stricmp( okmol[i].aT2, tAtom4 ) == 0 && stricmp( okmol[i].aT3, tAtom2 ) == 0 && stricmp( okmol[i].aT4, tAtom3 ) == 0 )
						  || ( stricmp( okmol[i].aT2, tAtom4 ) == 0 && stricmp( okmol[i].aT3, tAtom3 ) == 0 && stricmp( okmol[i].aT4, tAtom2 ) == 0 ) 
						  || ( stricmp( tAtom2, "X" ) == 0 && stricmp( okmol[i].aT3, tAtom3 ) == 0 && stricmp( okmol[i].aT4, tAtom4 ) == 0 )
						  || ( stricmp( tAtom2, "X" ) == 0 && stricmp( okmol[i].aT3, tAtom4 ) == 0 && stricmp( okmol[i].aT4, tAtom3 ) == 0 ) )
					  {
						  okmol[i].vn2 = tVn2;
					      okmol[i].gamma = tGamma;
						  okmol[i].n = tN;
					  }
				  }
			   }

	      }
        // end ( read improper dihedral block )

        // read non-bonding block
        if( IsNbBlock )
	   {
		  char tAtom[16];
		  strcpy( tAtom, "" ); // initialize

		  double tR = 0.0;
		  double tEps = 0.0;

		  char* pchar;

		  pchar = strtok( buffer, "-, "); // call first token ( separated by a hyphen or blank ) (atom type)
		  if( pchar != NULL )
	        { strcpy( tAtom, pchar );}

		  pchar = strtok( NULL, "-, "); // call second token ( Van der Waals radius )
	      if( pchar != NULL)
		     { tR = atof(pchar);} // Van der Waals radius in angstrom 

		  pchar = strtok( NULL, "-, "); // call second token ( Van der Waals well depth )
	      if( pchar != NULL)
		      { tEps = atof(pchar);} // Van der Waals well depth in kcal/mol

		  for( i = 0; i < nkmol.size() ; i++ )
		  {
			  if ( stricmp( nkmol[i].aT, tAtom ) ==0 )
			  {
			    	nkmol[i].R = tR;
			     	nkmol[i].eps = tEps;
			  }
		  }

	    }
	   // end (read non-bonding block)
      
	 } // end of while(1) (unfinite loop)

	      s.close(); // close a file stream

		  ts << tr("*>>>>>>>   General AMBER FF (GAFF) written in CHARMM FF style  <<<<<<<\n");
		  ts << tr("*>>>>>>>   for ") << molBaseName;
		  ts << tr(", generated by Avogadro                          <<<<<<<\n\n");
		  ts << "BONDS\n";
		  ts << "!\n";
          ts << "!V(bond) = Kb(b - b0)**2\n";
          ts << "!\n";
          ts << "!Kb: kcal/mole/A**2\n";
          ts << "!b0: A\n";
          ts << "!\n";
          ts << "!atom type Kb          b0\n";
          ts << "!\n";

		  for( i = 0; i < bkmol.size(); i++ )
	     {

		  ts << bkmol[i].aT1 << " "
             <<  bkmol[i].aT2 << "    "
			 << qSetRealNumberPrecision(9)
			 <<  bkmol[i].Kb << " "
			 <<  bkmol[i].req << endl;
   	     }

		  ts << endl; // insert empty line

		 // check if there are non-listed angles
	     vector<angleKind>  noListAk;

	     noListAk.clear();

	     for( i = 0; i < akmol.size(); i++ )
	     {
		     if( akmol[i].Kth == 0.0 && akmol[i].theq == 0.0 )
			    noListAk.push_back( akmol[i] );
	     }

	     vector<angleKind> akmol2;

	     if( !noListAk.empty() )
	     {
		    acquireWholeAngles( mol, &akmol2 );

		    for ( i = 0 ; i < noListAk.size() ; i++ )
		   {
			  int count = 0;

		  	  double tempR1 = 0.0;
			  double tempR2 = 0.0;
			  double tempTheta = 0.0;

			  for( int j = 0; j < akmol2.size() ; j++ )
			  {
				 if( ( strcmp( noListAk[i].aT1, akmol2[j].aT1 ) == 0 && strcmp( noListAk[i].aT2, akmol2[j].aT2 ) == 0 && strcmp( noListAk[i].aT3, akmol2[j].aT3 ) == 0 ) 
					|| ( strcmp( noListAk[i].aT1, akmol2[j].aT3 ) == 0 && strcmp( noListAk[i].aT2, akmol2[j].aT2 ) == 0 && strcmp( noListAk[i].aT3, akmol2[j].aT1 ) == 0 )  )
				  {
					tempR1 += akmol2[j].r1c;
					tempR2 += akmol2[j].r2c;
					tempTheta += akmol2[j].thc;

					count++;
				  }
			  }

			  if( count == 0 ) return -1;

			  noListAk[i].r1c = tempR1/(double)count;
			  noListAk[i].r2c = tempR2/(double)count;
			  noListAk[i].thc = tempTheta/(double)count;

			  /*
             // Calulate bond angle force constant for averaged distances and angles
             // of non-listed triple atom kinds in input molecule.
             // See J.Wang et al, J.Comp.Chem., 25, 1157-1174(2004)
             */
		      double D = pow( noListAk[i].r1c - noListAk[i].r2c, 2) / pow( noListAk[i].r1c + noListAk[i].r2c, 2);
		      double numerK = 143.9 * gaffAngleZ( noListAk[i].aN1 ) * gaffAngleC( noListAk[i].aN2 ) * gaffAngleZ( noListAk[i].aN3) * exp(-2*D);
		      double denomK = ( noListAk[i].r1c + noListAk[i].r2c ) * pow( noListAk[i].thc, 2);

			  noListAk[i].Kth = numerK/denomK;
			  noListAk[i].theq = noListAk[i].thc;
		     }

           for( i = 0 ; i < akmol.size() ; i ++ )
		   {
		      if( akmol[i].Kth == 0.0 && akmol[i].theq == 0.0 )
		      {
                   for( int j = 0 ; j < noListAk.size() ; j++ )
				   {
					 if( ( strcmp( akmol[i].aT1, noListAk[j].aT1) == 0 && strcmp( akmol[i].aT2, noListAk[j].aT2) == 0 && strcmp( akmol[i].aT3, noListAk[j].aT3) == 0 )
						 || ( strcmp( akmol[i].aT1, noListAk[j].aT3) == 0 && strcmp( akmol[i].aT2, noListAk[j].aT2) == 0 && strcmp( akmol[i].aT3, noListAk[j].aT2) == 0 ) )
					   {
						 akmol[i].Kth = noListAk[j].Kth;
					     akmol[i].theq = noListAk[j].theq;
					   }
			       }
		       }
		   }

		 } // end (check of non-listed angles)

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

		 for( i = 0; i < akmol.size(); i++ )
	     {
		    ts <<  akmol[i].aT1 << " "
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

	   for( i = 0; i < tkmol.size(); i++ )
	   {
		  ts <<  tkmol[i].aT1 << " "
             <<  tkmol[i].aT2 << " "
			 <<  tkmol[i].aT3 << " "
             <<  tkmol[i].aT4 << "    "
			 <<  tkmol[i].vn2 << " "
			 <<  tkmol[i].n << "  "
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
	   for( i = 0; i < okmol.size(); i++ )
	  {
		ts <<  okmol[i].aT1 << " "
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
	for( i = 0; i< nkmol.size() ; i++ )
	{

		ts << nkmol[i].aT << "     "
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

}