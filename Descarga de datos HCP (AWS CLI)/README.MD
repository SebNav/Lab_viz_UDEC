# 📥 Descarga de Datos del Human Connectome Project (HCP) usando AWS CLI

Este repositorio proporciona una guía paso a paso para instalar y configurar la interfaz de línea de comandos de AWS (AWS CLI) en Ubuntu 22.04 con el objetivo de descargar datos del **Human Connectome Project (HCP)** alojados en Amazon S3.

---

## 🔧 Requisitos

- Ubuntu 22.04
- Conexión a internet
- Cuenta en [ConnectomeDB](https://db.humanconnectome.org/)
- Credenciales de acceso S3 proporcionadas por HCP (AWS Access Key y Secret)

---

## ⚙️ Instalación de AWS CLI

### 1. Descargar e instalar AWS CLI v2

```bash
sudo apt update
sudo apt install unzip curl -y

curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install```bash


### 2. Verificar la instalación
