; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "CCPi (ImageJ Plugin)"
#define MyAppVersion "0.1"
#define MyAppPublisher "STFC"
#define MyAppURL "http://ccpforge.cse.rl.ac.uk/gf/project/iqa/"

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{033B5EFE-8F45-451C-81A1-CCCD32B783F2}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
;AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={pf}\
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
DisableDirPage=yes
OutputBaseFilename=CCPi ImageJ Plugin(v{#MyAppVersion})
Compression=lzma
SolidCompression=yes
WizardImageFile=CCPi_Logo.bmp
WizardImageStretch=no
WizardSmallImageFile=CCPi_Small_Logo.bmp

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
Source: "ImageJ\*"; DestDir: "{code:GetImageJDirectory}"; Flags: ignoreversion recursesubdirs createallsubdirs

; NOTE: Don't use "Flags: ignoreversion" on any shared system files


[Code]
var 
  ImageJVersionPage: TInputOptionWizardPage;
  ImageJInstallationDirPage: TInputDirWizardPage;

  procedure InitializeWizard;
  begin

   ImageJInstallationDirPage := CreateInputDirPage(wpWelcome, 'Select ImageJ Installation directory','',  'To continue, click Next. If you would like to select a different folder, click Browse.', False,'');
   ImageJInstallationDirPage.Add('');

   ImageJInstallationDirPage.Values[0] := GetPreviousData('ImageJDirectory', '');

  end;

  procedure RegisterPreviousData(PreviousDataKey: Integer);
begin
  { Store the settings so we can restore them next time }
  SetPreviousData(PreviousDataKey, 'ImageJDirectory', ImageJInstallationDirPage.Values[0]);
end;

function NextButtonClick(CurPageID: Integer): Boolean;
var
  I: Integer;
begin
  { Validate certain pages before allowing the user to proceed }
  if CurPageID = ImageJInstallationDirPage.ID then begin
    if ImageJInstallationDirPage.Values[0] = '' then
      begin
        MsgBox('You must a valid ImageJ Root Directory.', mbError, MB_OK);
        Result := False;
      end
    else if not FileExists( ImageJInstallationDirPage.Values[0]+'\ImageJ.exe') then
      begin
        MsgBox('You must a valid ImageJ Root Directory.', mbError, MB_OK);
        Result := False;
      end
    else
      Result := True;
  end else
    Result := True;
end;

function GetImageJDirectory(Param: String): String;
begin
  { Return the selected Installation Directory }
  Result := ImageJInstallationDirPage.Values[0];
end;


