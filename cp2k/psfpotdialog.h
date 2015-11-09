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

#ifndef PSFPOTDIALOG_H
#define PSFPOTDIALOG_H


#include <avogadro/molecule.h>

#include <QtCore/QSettings>
#include <QtGui/QDialog>
#include <QFileInfo>

#include "ui_psfpotdialog.h"

using namespace OpenBabel;

namespace Avogadro
{
   // Three struct used for building charmm-style parameter file
   struct bondKind
  {
	char    aT1[3];    // atom Type 1, two roman letters and '\0'
	char    aT2[3];	   // atom Type 2, two roman letters and '\0'
    double  Kb;     // force constant (kcal/mol). it is 0.0 initially
	double  req;    // equibrium distance (A).

  };

   struct angleKind
  {
	char  aT1[3]; // atom Type 1
	char  aT2[3]; // atom Type 2
	char  aT3[3]; // atom Type 3
	double  Kth;  // constant of bending motion (kcal/(mol radian^2))
	double  theq; // equibrium angle (degree)
	unsigned int aN1; // atomic number of 1. added in version 3
	unsigned int aN2;
	unsigned int aN3;
	double  r1c; // currnet bond length (A) between 1-2 atoms. added in version 3
	double  r2c; // currnet bond length (A) between 2-3 atoms. added in version 3
    double  thc; // current angle (degree). added in version 3
   };

   struct torKind
  {
	char  aT1[3];
	char  aT2[3];
	char  aT3[3];
	char  aT4[3];
	int   m;        // number of bond paths. ignoired.
	double   vn2;   // equals to Kchi for charmm pot (kcal/mol)
	double   gamma; // equals to delta for charmm pot (degree)
    double  n;      // the periodicity of torsion
   };

   struct oopKind
  {
	char  aT1[3];
	char  aT2[3];
	char  aT3[3];
	char  aT4[3];
	double   vn2;   // equals to Kchi for charmm pot (kcal/mol)
	double   gamma; // equals to delta for charmm pot (degree)
    double  n;      // the periodicity of torsion
   };

   struct nbKind
  {
	char aT[3];
	double  R;
	double eps;
   };

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

    private:
      Ui::PsfPotDialog ui;

	  QString psfPreviewPane();
	  QString potPreviewPane();

	 QString psfGaff();
	 QString potGaff();

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

	 void acquireWholeAngles( OBMol mol, std::vector<angleKind>* ak );
	 double gaffAngleC( unsigned int atomicNum );
	 double gaffAngleZ( unsigned int atomicNum );

	 bool Updated;

   public Q_SLOTS:
     void updatePreviewText();

   private Q_SLOTS:
     //! Button Slots
     void psfGenerateClicked();
	 void potGenerateClicked();

   protected:
	  QString saveInputFile(QString inputDeck, QString fileType, QString ext);
	  Molecule* m_molecule;
	  QString m_savePath;

  };

 }

#endif
