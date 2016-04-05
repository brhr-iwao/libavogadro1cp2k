/**********************************************************************
  Cp2kOutputDialog - Analyze CP2k Output

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

#ifndef CP2KOUTPUTDIALOG_H
#define CP2KOUTPUTDIALOG_H

#include <avogadro/molecule.h>
#include <avogadro/atom.h>
#include <avogadro/animation.h>

#include <openbabel/mol.h>
//#include <openbabel/generic.h> 

#include <QDialog>
#include <QTableWidget>
#include <QHeaderView>
//#include <QLineEdit>
//#include <QDoubleValidator>
#include <QTableWidgetSelectionRange>
#include <QPushButton>
#include <QLineEdit>
#include <QSlider>
#include <QCheckBox>
//#include <QGroupBox>
#include <QDebug>

#include <QFileInfo>
#include <QDir>
#include <QFileDialog>
//#include <QMessageBox>
#include <QIODevice>

#include "ui_cp2koutputdialog.h"
#include "ui_spectraplotdialog.h"

// Forward declaration of Avogadro::Molecule
namespace Avogadro {
  class Molecule;
}

using namespace OpenBabel;

namespace Avogadro
{
	class Cp2kOutputDialog : public QDialog
	{
	  Q_OBJECT

	   public:
          // Constructor
          explicit Cp2kOutputDialog( QWidget * parent = 0, Qt::WindowFlags flags = 0 );

          // Deconstructor
          virtual ~Cp2kOutputDialog();

		  void setMolecule(Molecule *molecule);

		public slots:

		private slots:
			void closeClicked();
			void loadFile();

			bool getVibrationWidgets();

			void createVibrationAnimation( );

        private:

		   Ui::Cp2kOutputDialog ui;

		   Molecule *m_molecule;

		   // vibration widgets of Avogadro
		   QWidget* vibrationDock;
		   QWidget* vibrationWidget;
		   QTableWidget* vibrationTable;
		   QLineEdit* editFilter;
		   QLabel* filterLabel;
		   QLabel* kmmolLabel;
		   QPushButton* spectraButton;
		   QSlider* scaleSlider;
		   QCheckBox* normalizeDispCheckBox;
		   QCheckBox* displayForcesCheckBox;
		   QCheckBox* animationSpeedCheckBox;
		   QPushButton* animationButton;
		   QPushButton* pauseButton;

		   // vibration data
		   std::vector< std::vector< vector3 > > vLx;
		   std::vector< double > vF;
		   std::vector< double > vI;

		   // vibration animation
		   Animation* m_animation;
		   std::vector<std::vector<Eigen::Vector3d> *> m_curFrames;
		   unsigned int m_framesPerStep;
		   qreal m_vibScale;
		   int iMode; // index of selected nomal mode

		   QString m_savePath;
           QString m_saveFilter;


  };// end class Cp2kOutputDialog


  class SpectraPlotDialog : public QDialog
  {
    Q_OBJECT

    public:
      explicit SpectraPlotDialog(QWidget *parent = 0, Qt::WindowFlags f = 0 );
      virtual ~SpectraPlotDialog();

    public slots:
      void setIRSpectra(Molecule *mol, std::vector< double > vF, std::vector< double > vI);

    private:
      // This member provides access to all ui elements
      Ui::SpectraPlotDialog ui;
      Molecule *m_molecule;

  };// end SpectraPlotDialog


} // end namespace Avogadro



#endif // #ifndef CP2KOUTPUTDIALOG_H