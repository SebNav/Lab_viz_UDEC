# 📥 Descarga de Datos del Human Connectome Project (HCP) usando AWS CLI

Este repositorio proporciona una guía paso a paso para instalar y configurar la interfaz de línea de comandos de AWS (AWS CLI) en Ubuntu 22.04 con el objetivo de descargar datos del **Human Connectome Project (HCP)** alojados en Amazon S3.

---

## 🔧 Requisitos

- Cuenta en [ConnectomeDB](https://db.humanconnectome.org/)
- **Credenciales de acceso S3 proporcionadas por HCP**  
  > 📌 Estas credenciales se obtienen al registrarte en el portal [ConnectomeDB](https://db.humanconnectome.org/), aceptar el acuerdo de uso de datos, y solicitar acceso a los datos restringidos. Una vez aprobado, recibirás un `Access Key ID` y un `Secret Access Key`.

---

## ⚙️ Instalación de AWS CLI

### 1. Descargar e instalar AWS CLI v2

```bash
sudo apt update
sudo apt install unzip curl -y

curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```


### 2. Verificar la instalación

```bash
aws --version
```

Deberías aparecer una salida similar a:


```bash
aws-cli/2.x.x Python/x.x.x Linux/x86_64

```


## 🔐 Configuración de credenciales

Una vez que tengas tus claves de acceso S3 proporcionadas por HCP (puedes solicitarlas registrándote en [ConnectomeDB](https://db.humanconnectome.org/), aceptando el acuerdo de uso de datos y solicitando acceso a los datos restringidos), configura el acceso con:

```bash
aws configure
```

se te pedira ingresar:

- AWS Access Key ID: tu clave pública (proporcionada por HCP)
- AWS Secret Access Key: tu clave secreta (proporcionada por HCP)
- Region: puedes escribir us-east-1 o dejar en blanco
- Output format: puedes dejar en blanco o usar json


## 📦 Descarga de datos HCP

Los datos del HCP están disponibles públicamente en el bucket s3://hcp-openaccess. Para descargar datos se debe colocar "aws s3" seguido de 'cp' o 'sync' dependiendo si se quiere descargar un archivo especifico o una carpeta completa, despues de esto se debe escribir la direccion del archivo o carpeta a descargar y por ultimo el lugar donde se guardara.

- Descargar un archivo específico

```bash
aws s3 cp s3://hcp-openaccess/HCP_1200/100307/T1w/T1w_acpc_dc_restore.nii.gz /home/prueba/T1w_acpc_dc_restore.nii.gz

```

- Descarga de un carpeta en especifico

```bash
aws s3 cp s3://hcp-openaccess/HCP_1200/100307/T1w/ /home/prueba/T1w/

```




