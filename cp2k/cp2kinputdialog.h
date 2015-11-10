/**********************************************************************
  Cp2kInputDialog - CP2k Input Dialog

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
	void readSettings(QSettings&);
    void writeSettings(QSettings&) const;

	void setModel(ConstraintsModel *model);

	// enum runType{ENERGY, ENERGY_FORCE, MD, GEO_OPT, MC};

  private:
    Ui::Cp2kInputDialog ui;
	
	QString m_projectName;
	// runType m_runType; 
	QString m_runType;
    // Molecule* m_molecule;
    ConstraintsModel *m_constraints;

	bool m_mmRadioChecked;
	bool m_qmRadioChecked;

	std::vector<QString> atomKind;

	QString generateInputDeck();

    void setAtomKind();


  public Q_SLOTS:
    void updatePreviewText();

    void mmRadioChecked();
	void qmRadioChecked();

  private Q_SLOTS:
    //! Button Slots
    void resetClicked();
    void generateClicked();

	void setProjectName();
	void setRunType(int);

  protected:
	QString saveInputFile(QString inputDeck, QString fileType, QString ext);

	Molecule* m_molecule;
    int m_multiplicity;
    int m_charge;
    QString m_savePath;

  };
}

#endif