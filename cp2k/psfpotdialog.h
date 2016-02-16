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

#ifndef PSFPOTDIALOG_H
#define PSFPOTDIALOG_H

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/residue.h>
#include <openbabel/forcefields/forcefieldgaff.h>
#include <openbabel/forcefields/forcefielduff.h>

#include <avogadro/molecule.h>


#include <QtCore/QSettings>
#include <QtGui/QDialog>
#include <QFileInfo>

#include "ui_psfpotdialog.h"

#define EPS 0.01 // for double-double comparison

using namespace std;
using namespace OpenBabel;

namespace Avogadro
{
   // Five structures used for building charmm-style parameter file
   struct bondKind
  {
	char    aT1[6];    // atom Type 1
	char    aT2[6];	   // atom Type 2
    double  Kb;     // force constant (kcal/mol). it is 0.0 initially
	double  req;    // equibrium distance (A).

  };

   struct angleKind
  {
	char  aT1[6]; // atom Type 1
	char  aT2[6]; // atom Type 2
	char  aT3[6]; // atom Type 3
	double  Kth;  // constant of bending motion (kcal/(mol radian^2))
	double  theq; // equibrium angle (degree)
   };

   struct torKind
  {
	char  aT1[6];
	char  aT2[6];
	char  aT3[6];
	char  aT4[6];
	int   m;        // number of bond paths. ignoired.
	double   vn2;   // equals to Kchi for charmm pot (kcal/mol)
	double   gamma; // equals to delta for charmm pot (degree) (phase)
    double  n;      // the periodicity of torsion (multiplicity)
   };

   struct oopKind
  {
	char  aT1[6];
	char  aT2[6];
	char  aT3[6];
	char  aT4[6];
	double   vn2;   // equals to Kchi for charmm pot (kcal/mol)
	double   gamma; // equals to delta for charmm pot (degree)
    double  n;      // the periodicity of torsion
   };

   struct nbKind
  {
	char aT[6];
	double  R;
	double eps;
  };

  
   // global variables to be stored with FF parameters
	static vector<bondKind>  bkmol;
	static vector<angleKind> akmol;
	static vector<torKind>   tkmol;
    static vector<oopKind>   okmol;
    static vector<nbKind>    nkmol;
   

  class ffpGaff : public OBForceFieldGaff
  {
	   public:
      
        bool SetupPointers()
		{
			   if( &( _bondcalculations ) == NULL) return false;

			   bkmol.clear();

			   vector<OBFFBondCalculationGaff>::iterator itb;

			   for (itb = _bondcalculations.begin(); itb != _bondcalculations.end(); ++itb) 
			   {
				   bondKind tempBk;

				   strcpy( tempBk.aT1, (*itb).a->GetType() );
				   strcpy( tempBk.aT2, (*itb).b->GetType() );
				   tempBk.Kb = (*itb).kr / KCAL_TO_KJ ;
				   tempBk.req = (*itb).r0;

				   bkmol.push_back( tempBk );

			   }


			   if( &( _anglecalculations ) == NULL ) return false;

			   akmol.clear();

			   vector<OBFFAngleCalculationGaff>::iterator ita;

			   for (ita = _anglecalculations.begin(); ita != _anglecalculations.end(); ++ita) 
			   {
				   angleKind tempAk;

				   strcpy( tempAk.aT1, (*ita).a->GetType() );
				   strcpy( tempAk.aT2, (*ita).b->GetType() );
				   strcpy( tempAk.aT3, (*ita).c->GetType() );
				   tempAk.Kth = (*ita).kth / KCAL_TO_KJ;
				   tempAk.theq = (*ita).theta0;

				   akmol.push_back( tempAk );
			   }

               if( &( _torsioncalculations ) == NULL ) return false;

			   tkmol.clear();

			   vector<OBFFTorsionCalculationGaff>::iterator itt;

			   for (itt = _torsioncalculations.begin(); itt != _torsioncalculations.end(); ++itt)
			   {
				   torKind tempTk;

				   strcpy( tempTk.aT1, (*itt).a->GetType() );
				   strcpy( tempTk.aT2, (*itt).b->GetType() );
				   strcpy( tempTk.aT3, (*itt).c->GetType() );
				   strcpy( tempTk.aT4, (*itt).d->GetType() );
				   tempTk.m = 0;
				   tempTk.vn2 = (*itt).vn_half / KCAL_TO_KJ;
				   tempTk.gamma = (*itt).gamma;
				   tempTk.n = (*itt).n;

				   tkmol.push_back( tempTk );
			   }


			   if( &( _oopcalculations ) == NULL ) return false;

			   okmol.clear();

			   vector<OBFFOOPCalculationGaff>::iterator ito;

			   for (ito = _oopcalculations.begin(); ito != _oopcalculations.end(); ++ito) 
			   {
				   oopKind tempOk;

				   strcpy( tempOk.aT1, (*ito).a->GetType() );
				   strcpy( tempOk.aT2, (*ito).b->GetType() );
				   strcpy( tempOk.aT3, (*ito).c->GetType() );
				   strcpy( tempOk.aT4, (*ito).d->GetType() );
				   tempOk.vn2 = (*ito).vn_half / KCAL_TO_KJ;
				   tempOk.gamma = (*ito).gamma;
				   tempOk.n = (*ito).n;

				   okmol.push_back( tempOk );
			   }

			   if( &( _ffvdwparams ) == NULL ) return false;

			   nkmol.clear();

			   for (int i = 0; i < _ffvdwparams.size() ; i++ )
			   {
				   nbKind tempNk;

				   strcpy( tempNk.aT, _ffvdwparams[i]._a.c_str() );

				  FOR_ATOMS_OF_MOL( a, _mol )
			      {
					 char* aType = a->GetType();

					 if( aType != NULL && strcmp( aType, tempNk.aT ) == 0)
					  {
						  tempNk.R = _ffvdwparams[i]._dpar[0];
				          tempNk.eps = _ffvdwparams[i]._dpar[1];

						  nkmol.push_back( tempNk );
					  }
			      }
		  
			   }

			   return true;

		    }


  }; 
  // end of class ffpGaff definition

  class ffpUff : OBForceFieldUFF
  {
	  public:
      
        bool SetupPointers()
		{
			   if( &( _bondcalculations ) == NULL) return false;

			   bkmol.clear();

			   vector<OBFFBondCalculationUFF>::iterator itb;

			   for (itb = _bondcalculations.begin(); itb != _bondcalculations.end(); ++itb ) 
			   {
				   bondKind tempBk;

				   strcpy( tempBk.aT1, (*itb).a->GetType() );
				   strcpy( tempBk.aT2, (*itb).b->GetType() );
				   tempBk.Kb = (*itb).kb / KCAL_TO_KJ ;
				   tempBk.req = (*itb).r0;

				   bkmol.push_back( tempBk );
			   }


			   if( &( _anglecalculations ) == NULL ) return false;

			   akmol.clear();

			   vector<OBFFAngleCalculationUFF>::iterator ita;

			   for (ita = _anglecalculations.begin(); ita != _anglecalculations.end(); ++ita ) 
			   {
				   angleKind tempAk;

				   strcpy( tempAk.aT1, (*ita).a->GetType() );
				   strcpy( tempAk.aT2, (*ita).b->GetType() );
				   strcpy( tempAk.aT3, (*ita).c->GetType() );
				   tempAk.Kth = (*ita).ka / KCAL_TO_KJ;
				   tempAk.theq = (*ita).theta0;

				   akmol.push_back( tempAk );
			   }

			   if( &( _torsioncalculations ) == NULL ) return false;

			   tkmol.clear();

			   vector<OBFFTorsionCalculationUFF>::iterator itt;

			   for (itt = _torsioncalculations.begin(); itt != _torsioncalculations.end(); ++itt )
			   {
				   torKind tempTk;

				   strcpy( tempTk.aT1, (*itt).a->GetType() );
				   strcpy( tempTk.aT2, (*itt).b->GetType() );
				   strcpy( tempTk.aT3, (*itt).c->GetType() );
				   strcpy( tempTk.aT4, (*itt).d->GetType() );
				   tempTk.m = 0;
				   tempTk.vn2 = (*itt).V / KCAL_TO_KJ;

				      // j,k are both sp3
					   if (((*itt).b->GetType())[2] == '3' && ((*itt).c->GetType())[2] == '3' ) 
					   {
						   tempTk.gamma = 180;
						   tempTk.n = 3;
					   }

					   else if ((((*itt).b->GetType())[2] == '2' && ((*itt).c->GetType())[2] == '3')
						   || (((*itt).b->GetType())[2] == '3' && ((*itt).c->GetType())[2] == '2')
						   || ( strcmp( (*itt).b->GetType(),"C_R") == 0 && ((*itt).c->GetType())[2] == '3')
						   || (((*itt).b->GetType())[2] == '3' && strcmp((*itt).c->GetType(), "C_R")== 0 ) )
					   {
						   tempTk.gamma = 0;
						   tempTk.n = 6;
					   }

					   // j,k are both sp2
					   else if ( (((*itt).b->GetType())[2] == '2' || strcmp((*itt).b->GetType(), "C_R") == 0 )
						   && ( ((*itt).c->GetType())[2] == '2' || strcmp((*itt).c->GetType(), "C_R") == 0  ) )
					   {
						   tempTk.gamma = 180;
						   tempTk.n = 2;
					   }

					   else
					   {
						   tempTk.gamma = 0;
						   tempTk.n = 1;
					   }


				   tkmol.push_back( tempTk );
			   }

			   if( &( _oopcalculations ) == NULL ) return false;

			   okmol.clear();

			   vector<OBFFOOPCalculationUFF>::iterator ito;

			   for (ito = _oopcalculations.begin(); ito != _oopcalculations.end(); ++ito ) 
			   {
				   oopKind tempOk;

				   strcpy( tempOk.aT1, (*ito).a->GetType() );
				   strcpy( tempOk.aT2, (*ito).b->GetType() );
				   strcpy( tempOk.aT3, (*ito).c->GetType() );
				   strcpy( tempOk.aT4, (*ito).d->GetType() );
				   tempOk.vn2 = (*ito).koop / KCAL_TO_KJ;
				   tempOk.gamma = 0;
				   tempOk.n = 1;

				   okmol.push_back( tempOk );
			   }

			  if( &(_ffparams) == NULL ) return false;

			  nkmol.clear();

			  for (int i = 0; i < _ffparams.size() ; i++ )
			  {
				   nbKind tempNk;

				   strcpy( tempNk.aT, _ffparams[i]._a.c_str() ); //Atom type

				  FOR_ATOMS_OF_MOL( a, _mol )
			      {
					 char* aType = a->GetType();

					 if( aType != NULL && strcmp( aType, tempNk.aT ) == 0)
					  {
						  tempNk.R = _ffparams[i]._dpar[2]; //x1
				          tempNk.eps = _ffparams[i]._dpar[3]; //D1

						  nkmol.push_back( tempNk );
					  }

			      }

			   }

	       return true;

		}

  };
  // end of class ffpUff definition

  class Molecule;
  class PsfPotDialog : public QDialog
  {
    Q_OBJECT

    public:
     explicit PsfPotDialog(QWidget *parent = 0, Qt::WindowFlags f = 0 );
     ~PsfPotDialog();

	 void setMolecule(Molecule *molecule);
	 void readSettings(QSettings&);
     void writeSettings(QSettings&) const;

	 QString psfGaff();
	 QString potGaff();

	 QString psfUff();
	 QString potUff();

    private:
      Ui::PsfPotDialog ui;

	  QString m_forceField;

	  QString psfPreviewPane();
	  QString potPreviewPane();

	 static int chrs2int( const char* chrs );

	 static bool bondTypeComp( const bondKind& left, const bondKind& right );
	 static bool bondTypeEqual(const bondKind& left, const bondKind& right );
	 static bool angleTypeComp( const angleKind& left, const angleKind& right );
	 static bool angleTypeEqual( const angleKind& left, const angleKind& right );
	 static bool torTypeComp( const torKind& left, const torKind& right );
	 static bool torTypeEqual( const torKind& left, const torKind& right );
	 static bool oopTypeComp( const oopKind& left, const oopKind& right );
	 static bool oopTypeEqual( const oopKind& left, const oopKind& right );
	 static bool nbTypeComp( const nbKind& left, const nbKind& right );
	 static bool nbTypeEqual( const nbKind& left, const nbKind& right );

	 bool Updated;

	 int connAtoms( OBAtom* a );

   public Q_SLOTS:
     void updatePreviewText();

   private Q_SLOTS:
     //! Button Slots
     void psfGenerateClicked();
	 void potGenerateClicked();
	 void closeClicked();

	 void setForceField(int n);

   protected:
	  QString saveInputFile(QString inputDeck, QString fileType, QString ext);
	  Molecule* m_molecule;
	  QString m_savePath;

  };

 }

#endif
