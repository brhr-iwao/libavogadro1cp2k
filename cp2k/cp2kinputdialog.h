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
#include "../constraintsmodel.h"


namespace Avogadro
{
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

	void setModel(ConstraintsModel *model);

  private:
    Ui::Cp2kInputDialog ui;
	
	// basic tab
	QString m_projectName;
	QString m_runType;

	// int m_viewAtomUid;
	bool m_viewAtomUid;
	// Qt::CheckState m_viewAtomUid;

    bool m_mmRadioChecked;
	bool m_qmRadioChecked;
	bool m_qmmmRadioChecked;

	// MM tab
	double m_emaxSpline;
	QString m_ewaldType;

	// QM tab
	QString m_qmMethod;

	// DFT tab
	QString m_basisSet;
	QString m_functional;

	int m_nMultiGrid;
	int m_cutOff;

    ConstraintsModel *m_constraints;

	std::vector<QString> atomKind;

	QString generateInputDeck();

    void setAtomKindMol();

	QString potentialName( QString atomType );


  public Q_SLOTS:
    void updatePreviewText();

    void mmRadioChecked();
	void qmRadioChecked();
	void qmmmRadioChecked();

	void setAtomLabelUid();

  private Q_SLOTS:
    //! Button Slots
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

	// DFT tab
    void setBasisSet(int);
	void setFunctional(int);
	void setNMultiGrid(int);
	void setCutOff(int);

  protected:
	QString saveInputFile(QString inputDeck, QString fileType, QString ext);

	Molecule* m_molecule;
    int m_multiplicity;
    int m_charge;
    QString m_savePath;

  };
}

#endif