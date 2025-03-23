pipeline{
    agent{
        docker{
            image 'debian:lastest'
            args '-u root:root'
        }
        environment{
            GIT_REPO = 'https://github.com/malaku'
            SSH_PASS = 'temppwd'
        }
    }
    stages{
        stage('Prepare Environment'){
            steps {
                sh 'apt-get update && apt-get install -y sshpass git'
            }
        }
            
        stage('Clone Repos'){
            steps{
                def repos = ['Vehicle-and-Driver-Model', 'FW1', 'FW2', 'RW1', 'RW2', 'WAC1', 'WAC2', 'WAC3', 'WAC4']
                for (repo in repos){
                    sh "git clone ${env.GIT_REPO}/${repo}.git"
                }
            }
        }
        stage('Run on Beaglebones'){
            steps {
                script {
                    def projects = [
                        [ board: 'beaglebone1', project: 'FW1'],
                        [ board: 'beaglebone2', project: 'FW2'],
                        [ board: 'beaglebone3', project: 'RW1'],
                        [ board: 'beaglebone4', project: 'RW2'],
                        [ board: 'beaglebone5', project: 'WAC1'],
                        [ board: 'beaglebone6', project: 'WAC2'],
                        [ board: 'beaglebone7', project: 'WAC3'],
                        [ board: 'beaglebone8', project: 'WAC4'],
                        [ board: 'beaglebone9', project: 'Vehicle-and-Driver-Model']
                    ]
                    for(entry in projects){
                        def remoteCmd = """
                            cd ${entry.project} && \
                            chmod +x init_can.sh build.sh run.sh && \
                            ./init_can.sh && \
                            ./build.sh && \
                            ./run.sh
                            """
                            sh "sshpass -p '${env.SSH_PASS}' ssh -o StrictHostKeyChecking=no debian@${entry.board} '${remoteCmd}'"
                    }
                }
            }
        }
        
        
}