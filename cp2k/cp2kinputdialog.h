/**********************************************************************
  Cp2kInputDialog - CP2k Input Dialog

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

#ifndef CP2KINPUTDIALOG_H
#define CP2KINPUTDIALOG_H


#include <avogadro/molecule.h>

#include <QtCore/QSettings>
#include <QtGui/QDialog>
#include <QFileInfo>

#include "ui_cp2kinputdialog.h"

namespace Avogadro
{

  struct selatom // selected atom
  {
	  QString element;
	  unsigned long uid; // atom unique id + 1
	  unsigned int idx;
  };

  struct linkatoms // qm-mm boundary atoms
  {
	  unsigned long mmuid; // unique id + 1 for mm atom
	  unsigned long qmuid;
      unsigned int  mmidx;
	  unsigned int  qmidx;
  };

  class Molecule;
  class Cp2kInputDialog : public QDialog
  {
  Q_OBJECT

  public:
    explicit Cp2kInputDialog(QWidget *parent = 0, Qt::WindowFlags f = 0 );
    ~Cp2kInputDialog();

	void setMolecule(Molecule *molecule);

    void writeSettings(QSettings&) const;
	void readSettings(QSettings&);

  private:
    // Variables
	// Common
    Ui::Cp2kInputDialog ui;

    std::vector<QString> atomKindMol;
	std::vector<selatom> selAtoms;
	std::vector<QString> selAtomsKind;
	std::vector<linkatoms> linkAtoms;
	
	// basic tab
	QString m_projectName;
	QString m_runType;

	bool m_viewAtomUid;

    bool m_mmRadioChecked;
	bool m_qmRadioChecked;
	bool m_qmmmRadioChecked;

	// MM tab
	double m_emaxSpline;
	QString m_ewaldType;

	// QM tab
	QString m_qmMethod;
	QString m_scfGuess;
	int m_maxSCF;

	// DFT tab
	QString m_basisSet;
	QString m_functional;

	int m_nMultiGrid;
	int m_cutOff;

	// DFTB tab
	bool m_isScc;
	bool m_dispersion;

	// SE tab
	QString m_seMethod;

	// QMMM tab
	QString m_qmmmElectrostaticCoupling;

	// MD tab
	QString m_ensemble;

   // Functions
	// Common
	QString generateInputDeck();
    void setAtomKindMol();
	QString potentialName( QString atomType );
	static bool qstringEqual(const QString& left, const QString& right);
	void setSelAtoms();
	QString unitCell();
	void setLinkAtoms();

  public Q_SLOTS:
    void updatePreviewText();

    void mmRadioChecked();
	void qmRadioChecked();
	void qmmmRadioChecked();

  private Q_SLOTS:
    // Button Slots
    void resetClicked();
    void generateClicked();
	void closeClicked();

	// basic tab
	void setProjectName();
	void setRunType(int);

	// MM tab
    void setEmaxSpline(double);
	void setEwaldType(int);

	// QM tab
	void setQmMethod(int);
	void setCharge(int);
	void setMultiplicity(int);
	void setMaxSCF(int);
	void setSCFGuess(int);

	// DFT tab
    void setBasisSet(int);
	void setFunctional(int);
	void setNMultiGrid(int);
	void setCutOff(int);

	// DFTB tab
	void setDftbSCCType(int);
	void setDftbDisp(int);

	// SE tab
	void setSEMethod(int);

	// QMMM tab
	void setQMMMElectrostaticCoupling(int);

	// MD tab
	void setEnsemble(int);

  protected:
	QString saveInputFile(QString inputDeck, QString fileType, QString ext);

	Molecule* m_molecule;
    int m_multiplicity;
    int m_charge;
    QString m_savePath;

  };
}

#endif